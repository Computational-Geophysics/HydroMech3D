#include "Mesh.h"                      // for Mesh and ElementData
// Make sure the header includes the full definition of BigWhamIO, not just a forward declaration.

#include <fstream>
#include <cmath>   // for fmax, etc.
#include <cassert>
#include <memory>
#include <iostream>
#include <array>

#include "ExportBackgroundStress.h"
#include "core/be_mesh.h"                     // for CreateMeshFromVect
#include "elements/rectangle.h"               // for Rectangle<0>
#include "elements/point.h"                   // for Point<3>
#include "core/elastic_properties.h"          // for ElasticProperties
#include "hmat/square_matrix_generator.h"     // for BieMatrixGenerator
#include "hmat/hierarchical_representation.h" // for HRepresentationRectangularMatrix
#include "hmat/hmatrix/hmat.h"                // for Hmat
#include "io/bigwham_io_helper.h"
#include "elasticity/bie_elastostatic.h"

namespace {
struct Vec3 { double x, y, z; };
static Vec3 cross(const Vec3 &a, const Vec3 &b) {
    return {a.y*b.z - a.z*b.y,
            a.z*b.x - a.x*b.z,
            a.x*b.y - a.y*b.x};
}
static Vec3 normal(const double *p0, const double *p1, const double *p2) {
    Vec3 e1{p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
    Vec3 e2{p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
    return cross(e1, e2);
}
static bool needsFlip(const Vec3 &n, int sign, char axis) {
    double axisVal = 0.0, offMag = 0.0;
    switch (axis) {
      case 'x':
        axisVal = n.x;
        offMag  = std::sqrt(n.y*n.y + n.z*n.z);
        break;
      case 'y':
        axisVal = n.y;
        offMag  = std::sqrt(n.x*n.x + n.z*n.z);
        break;
      case 'z':
        axisVal = n.z;
        offMag  = std::sqrt(n.x*n.x + n.y*n.y);
        break;
    }
    double axisMag = std::abs(axisVal);
    if (offMag > 1e-8 * axisMag)
        std::cerr << "[ExportBackgroundStress] Warning: "
                  << axis << "-panel normal not axis aligned\n";
    return sign * axisVal < 0;
}
}

namespace EQSim{
struct ExportBackgroundStress::Impl {
    int n_openMP_threads_;
    std::string kernel_name_;
    int bs_=3;

    std::unique_ptr<BigWhamIO> hmat_io_;

    std::vector<std::unique_ptr<il::Array2D<double>>> fullRankBlocks_;
    std::vector<std::unique_ptr<bigwham::LowRank<double>>> lowRankBlocks_;
    struct FR { int i0,i1,j0,j1; };
    struct LR { int i0,i1,j0,j1; };
    std::vector<FR> frBlocks_;
    std::vector<LR> lrBlocks_;
    std::vector<int> perm0_, perm1_;
    std::vector<double> rec_coords_;

    std::shared_ptr<bigwham::HRepresentation> hr_;
    std::shared_ptr<bigwham::Mesh> mesh_src_;
    std::shared_ptr<bigwham::Mesh> mesh_rec_;
    std::shared_ptr<bigwham::Mesh> mesh_; // for square matrix
    std::shared_ptr<bigwham::Mesh> mesh_obs_;

    std::shared_ptr<bigwham::BieKernel<double>> ker_obj_; // BieKernel description for the matrix operator
    std::shared_ptr<bigwham::BieKernel<double>> ker_obs_u_; // BieKernel description for computing observation of potential/displacement (u) at points
    std::shared_ptr<bigwham::BieKernel<double>> ker_obs_q_; // BieKernel description for computing observation of 'flux'/'stress' (q) at points

    Impl(const Mesh &faultMesh, double E, double nu)
{
    // 1. Determine observation grid dimensions and layer index
    // (Assume we can get NX, NY, NZ from the simulation's solid matrix or global properties)

    // 2. Prepare source (fault) mesh data for BigWham
    arma::mat coordsMat = faultMesh.getCoordinates();            // [N_nodes x 3] matrix
    arma::imat connMat  = faultMesh.getNodesConnettivity();      // [N_elts x m] (m=4 for rect, 3 for tri)
    arma::uword Nelts = faultMesh.getNumberOfElts();
    // Determine element type from connectivity
    int nodesPerElt = connMat.n_cols; 
    std::string kernelName;
    if (nodesPerElt == 4) {
        kernelName = "3DR0-3DR0-H";   // Rectangle elements:contentReference[oaicite:16]{index=16}
    } else if (nodesPerElt == 3) {
        kernelName = "3DT0-3DT0-H";   // Triangle elements
    } else {
        throw std::runtime_error("Unsupported element type in fault mesh");
    }
    // Populate coordinate and connectivity vectors for BigWham (source)
    std::vector<double> coor_src;
    coor_src.reserve(coordsMat.n_rows * 3);
    for (arma::uword i = 0; i < coordsMat.n_rows; ++i) {
        coor_src.push_back(coordsMat(i,0));
        coor_src.push_back(coordsMat(i,1));
        coor_src.push_back(coordsMat(i,2));
    }
    std::vector<int> conn_src;
    conn_src.reserve(connMat.n_rows * nodesPerElt);
    for (arma::uword e = 0; e < connMat.n_rows; ++e) {
        for (int n = 0; n < nodesPerElt; ++n) {
            conn_src.push_back(connMat(e,n));
        }
    }

    // 3. Construct receiver mesh for vertical snapshots at X = –10, 0, and 10 m
    std::vector<double> xPlanes = { 0.0 };
    int NY_ = 81;
    int NZ_ = 81;
    double y_min = -18.0, y_max = 18.0;
    double z_min = -18.0, z_max = 18.0;
    double dy = (y_max - y_min) / (NY_ - 1);
    double dz = (z_max - z_min) / (NZ_ - 1);
    double eps = 0.25 * std::min(dy, dz);
    int recNodesPerElt = nodesPerElt;  // 4 for rectangles, 3 for triangles

    // Buffers for BigWham mesh
    std::vector<double> coor_rec;
    std::vector<int>    conn_rec;
    // Buffer for exactly one (x,y,z) centroid per panel
    std::vector<double> centroids;

    coor_rec.reserve(xPlanes.size() * NY_ * NZ_ * recNodesPerElt * 3);
    conn_rec.reserve(xPlanes.size() * NY_ * NZ_ * recNodesPerElt);
    centroids.reserve(xPlanes.size() * NY_ * NZ_ * 3);

    int nodeIndexOffset = 0;
    for (double X : xPlanes) {
      for (int j = 0; j < NY_; ++j) {
        double Y = y_min + j * dy;
        for (int k = 0; k < NZ_; ++k) {
          double Z = z_min + k * dz;

          //
          // -- X-normal panel (plane ⟂ X-axis, in Y–Z) --
          //
          // Record its centroid
          centroids.push_back(X);
          centroids.push_back(Y);
          centroids.push_back(Z);

          if (nodesPerElt == 4) {
            // rectangle: 4 corners
            std::array<double,12> verts = {
              X,      Y+eps, Z+eps,
              X,      Y-eps, Z+eps,
              X,      Y-eps, Z-eps,
              X,      Y+eps, Z-eps
            };
            Vec3 n = normal(&verts[0], &verts[3], &verts[6]);
            coor_rec.insert(coor_rec.end(), verts.begin(), verts.end());
            int signWanted = (X < 0.0 ? -1 : 1);
            if (needsFlip(n, signWanted, 'x')) {
              conn_rec.insert(conn_rec.end(), {
                nodeIndexOffset+0,
                nodeIndexOffset+3,
                nodeIndexOffset+2,
                nodeIndexOffset+1
              });
            } else {
              conn_rec.insert(conn_rec.end(), {
                nodeIndexOffset+0,
                nodeIndexOffset+1,
                nodeIndexOffset+2,
                nodeIndexOffset+3
              });
            }
            nodeIndexOffset += 4;
          } else {
            // triangle: 3 corners
            std::array<double,9> verts = {
              X,      Y+eps, Z+eps,
              X,      Y-eps, Z+eps,
              X,      Y+eps, Z-eps
            };
            Vec3 n = normal(&verts[0], &verts[3], &verts[6]);
            coor_rec.insert(coor_rec.end(), verts.begin(), verts.end());
            int signWanted = (X < 0.0 ? -1 : 1);
            if (needsFlip(n, signWanted, 'x')) {
              conn_rec.insert(conn_rec.end(), {
                nodeIndexOffset+0,
                nodeIndexOffset+2,
                nodeIndexOffset+1
              });
            } else {
              conn_rec.insert(conn_rec.end(), {
                nodeIndexOffset+0,
                nodeIndexOffset+1,
                nodeIndexOffset+2
              });
            }
            nodeIndexOffset += 3;
          }
          

          //
          // -- Y-normal panel (plane ⟂ Y-axis, in X–Z) --
          //
          // Record its centroid (must match exactly one panel per centroid)
          centroids.push_back(X);
          centroids.push_back(Y);
          centroids.push_back(Z);

          if (nodesPerElt == 4) {
            // rectangle: 4 corners
            std::array<double,12> verts = {
              X+eps,  Y, Z+eps,
              X-eps,  Y, Z+eps,
              X-eps,  Y, Z-eps,
              X+eps,  Y, Z-eps
            };
            Vec3 n = normal(&verts[0], &verts[3], &verts[6]);
            coor_rec.insert(coor_rec.end(), verts.begin(), verts.end());
            int signWanted = (Y > 0.0 ? 1 : -1);
            if (needsFlip(n, signWanted, 'y')) {
              conn_rec.insert(conn_rec.end(), {
                nodeIndexOffset+0,
                nodeIndexOffset+3,
                nodeIndexOffset+2,
                nodeIndexOffset+1
              });
            } else {
              conn_rec.insert(conn_rec.end(), {
                nodeIndexOffset+0,
                nodeIndexOffset+1,
                nodeIndexOffset+2,
                nodeIndexOffset+3
              });
            }

            nodeIndexOffset += 4;
          }
          else {
            // triangle panels
            std::array<double,9> verts = {
              X+eps, Y, Z+eps,
              X-eps, Y, Z+eps,
              X+eps, Y, Z-eps
            };
            Vec3 n = normal(&verts[0], &verts[3], &verts[6]);
            coor_rec.insert(coor_rec.end(), verts.begin(), verts.end());
            int signWanted = (Y > 0.0 ? 1 : -1);
            if (needsFlip(n, signWanted, 'y')) {
              conn_rec.insert(conn_rec.end(), {
                nodeIndexOffset+0,
                nodeIndexOffset+2,
                nodeIndexOffset+1
              });
            } else {
              conn_rec.insert(conn_rec.end(), {
                nodeIndexOffset+0,
                nodeIndexOffset+1,
                nodeIndexOffset+2
              });
            }
            nodeIndexOffset += 3;
          }

          //
          // -- Z-normal panel (plane ⟂ Z-axis, in X–Y) --
          //
          centroids.push_back(X);
          centroids.push_back(Y);
          centroids.push_back(Z);

          if (nodesPerElt == 4) {
            std::array<double,12> verts = {
              X+eps, Y+eps, Z,
              X-eps, Y+eps, Z,
              X-eps, Y-eps, Z,
              X+eps, Y-eps, Z
            };
            Vec3 n = normal(&verts[0], &verts[3], &verts[6]);
            coor_rec.insert(coor_rec.end(), verts.begin(), verts.end());
            int signWanted = (Z < 0.0 ? -1 : 1);
            if (needsFlip(n, signWanted, 'z')) {
              conn_rec.insert(conn_rec.end(), {
                nodeIndexOffset+0,
                nodeIndexOffset+3,
                nodeIndexOffset+2,
                nodeIndexOffset+1
              });
            } else {
              conn_rec.insert(conn_rec.end(), {
                nodeIndexOffset+0,
                nodeIndexOffset+1,
                nodeIndexOffset+2,
                nodeIndexOffset+3
              });
            }
            nodeIndexOffset += 4;
          } else {
            std::array<double,9> verts = {
              X+eps, Y+eps, Z,
              X-eps, Y+eps, Z,
              X+eps, Y-eps, Z
            };
            Vec3 n = normal(&verts[0], &verts[3], &verts[6]);
            coor_rec.insert(coor_rec.end(), verts.begin(), verts.end());
            int signWanted = (Z < 0.0 ? -1 : 1);
            if (needsFlip(n, signWanted, 'z')) {
              conn_rec.insert(conn_rec.end(), {
                nodeIndexOffset+0,
                nodeIndexOffset+2,
                nodeIndexOffset+1
              });
            } else {
              conn_rec.insert(conn_rec.end(), {
                nodeIndexOffset+0,
                nodeIndexOffset+1,
                nodeIndexOffset+2
              });
            }
            nodeIndexOffset += 3;
          }

        }
      }
    }

    // 4. Build the BigWham receiver mesh from coor_rec/conn_rec
    // 4. Initialize BigWhamIO to build the H-matrix for slip->traction mapping
    std::vector<double> properties = {E, nu};
    int nThreads = 8;  // use 8 threads (could also use omp_get_max_threads())
    std::cout << "Receiver mesh: " 
          << (coor_rec.size()/3) << " nodes, "
          << (conn_rec.size()/nodesPerElt) << " elements\n";
    assert(!coor_rec.empty() && !conn_rec.empty());
    std::cout << "Source mesh: " 
          << (coor_src.size()/3) << " nodes, "
          << (conn_src.size()/4) << " elements\n";
    assert(coor_src.size() % 3 == 0);
    assert(conn_src.size() % 4 == 0);
    std::cout << "Initializing BigWhamIO with " << nThreads << " threads..." << std::endl;
    //hmat_io_ = std::make_unique<BigWhamIO>(coor_src, conn_src, coor_rec, conn_rec, kernelName, properties, nThreads);
    bigwham::ElasticProperties elas(properties[0], properties[1]);
    il::int_t dof_dimension_ = 3;
    il::int_t spatial_dimension_ = 3;
    il::int_t flux_dimension_ = 6 ;
    using src_elem = bigwham::Rectangle<0>;
    using rec_elem = bigwham::Rectangle<0>;
    mesh_src_ = bigwham::CreateMeshFromVect<src_elem>(
            spatial_dimension_, /* num vertices */ 4, coor_src, conn_src);
    mesh_rec_ = bigwham::CreateMeshFromVect<rec_elem>(
    spatial_dimension_,/*num vertices per element*/ nodesPerElt,
    coor_rec, conn_rec);

    ker_obj_ = std::make_shared<bigwham::BieElastostatic<src_elem, rec_elem, bigwham::ElasticKernelType::T>>(
            elas, spatial_dimension_);
    mesh_ = mesh_src_;
    rec_coords_ = std::move(centroids);//rec_coords_ = coor_rec;

    // observations....
    using ObsType = bigwham::Point<3>;
    //ker_obs_u_=std::make_shared<bigwham::BieElastostatic<src_elem, ObsType, bigwham::ElasticKernelType::T>>(
    //       elas, spatial_dimension_);
    //ker_obs_q_=std::make_shared<bigwham::BieElastostatic<src_elem,ObsType, bigwham::ElasticKernelType::W>>(
    //        elas, spatial_dimension_);
    int leafSize = 100;
    double eta     = 3;
    //hmat_io_->BuildPattern(leafSize, eta);
    std::cout << "--------------------\n";
    std::cout << "Hierarchical representation creation ...\n";
    hr_ = HRepresentationRectangularMatrix(mesh_src_, mesh_rec_, leafSize, eta);
    std::cout << "Hierarchical representation complete.\n";
    // **2)** Populate the H‐matrix with ACA compression
    double epsACA  = 1e-4;
    //hmat_io_->BuildHierarchicalMatrix(leafSize, eta, epsACA);
    std::cout << "--------------------\n";
    std::cout << "Populating Hierarchical matrix ...\n";
    bigwham::BieMatrixGenerator<double> M(mesh_src_, mesh_rec_, ker_obj_, hr_);

    //hmat_ = std::make_shared<bigwham::Hmat<double>>(M, epsACA, n_openMP_threads_);
    auto &patt = hr_->pattern_;
    bs_ = 3;  // e.g. 3 for 3DR0-3DR0-H
    perm0_.assign(hr_->permutation_0_.begin(), hr_->permutation_0_.end());
    perm1_.assign(hr_->permutation_1_.begin(), hr_->permutation_1_.end());

    // --- Manual low-rank block build (no OpenMP, with debug) ---
    lowRankBlocks_.clear();
    lrBlocks_.clear();
    lowRankBlocks_.reserve(patt.n_LRB);
    lrBlocks_.reserve(patt.n_LRB);

    for (int k = 0; k < static_cast<int>(patt.n_LRB); ++k) {
        // extract the row/col ranges from the pattern
        int i0 =        patt.LRB_pattern(1, k);
        int i1 =        patt.LRB_pattern(3, k);
        int j0 =        patt.LRB_pattern(2, k);
        int j1 =        patt.LRB_pattern(4, k);

        std::cerr << "[ExportBackgroundStress] Building LR block #" 
                  << k << "  i:[" << i0 << "," << i1 
                  << "]  j:[" << j0 << "," << j1 << "]\n";

        try {
            // exactly as Hmatrix_bigwham does:
            auto lra = bigwham::adaptiveCrossApproximation<3>(
                         M,
                         {i0, i1},
                         {j0, j1},
                         epsACA
                       );
            // store it
            
            lowRankBlocks_.emplace_back(std::move(lra));
            lrBlocks_.push_back({i0, i1, j0, j1});
        }
        catch (const std::exception &ex) {
            std::cerr << "[ExportBackgroundStress] *** ACA failed at block #" 
                      << k << "  exception: " << ex.what() << "\n"
                      << "         i0="<<i0<<", i1="<<i1
                      <<", j0="<<j0<<", j1="<<j1<<"\n";
            throw;  // re-throw so you still see the crash
        }
    }
    std::cout << "Built " << lowRankBlocks_.size() 
              << " low-rank blocks.\n";

    // 4) full-rank blocks via direct set
    fullRankBlocks_.clear();
    frBlocks_.clear();
    fullRankBlocks_.reserve(patt.n_FRB);
    frBlocks_.reserve(patt.n_FRB);
    for (int k = 0; k < patt.n_FRB; ++k) {
        int i0 = patt.FRB_pattern(1, k);
        int i1 = patt.FRB_pattern(3, k);
        int j0 = patt.FRB_pattern(2, k);
        int j1 = patt.FRB_pattern(4, k);
        int ni = bs_ * (i1 - i0);
        int nj = bs_ * (j1 - j0);
        std::cerr << "[ExportBackgroundStress] Building FR block #" 
                  << k << "  i:[" << i0 << "," << i1 
                  << "]  j:[" << j0 << "," << j1 << "]\n";
        auto A = std::make_unique<il::Array2D<double>>(ni, nj);
        M.set(i0, j0, il::io, A->Edit());
        {
          fullRankBlocks_.emplace_back(std::move(A));
          frBlocks_.push_back({i0,i1,j0,j1});
        }
    }
    std::cout << "Built " << fullRankBlocks_.size() 
              << " full-rank blocks.\n";
    // 5) mark as built
    std::cout << "Hierarchical matrix population complete.\n";
    std::cout << "--------------------\n";
    std::cout << "BigWhamIO initialized successfully." << std::endl;
}



void ComputeAndExport(
    const std::vector<double> &slip,
    const std::string        &filename)
{
    // 1) copy slip into arma::vec
    arma::vec slip_vec(slip.size());
    for (size_t i = 0; i < slip.size(); ++i)
        slip_vec(i) = slip[i];

    // 2) apply H-matrix mat-vec
    arma::vec traction = multiply(slip_vec);

    // 3) open file and write header (with newline!)
    std::ofstream out(filename);
    out << "x,y,z,tx,ty,tz\n";

    // 4) loop over centroids
    int bs   = bs_;                          // DOFs per collocation
    int Npts = static_cast<int>(rec_coords_.size() / 3);
    for (int e = 0; e < Npts; ++e) {
        // fetch centroid
        double x  = rec_coords_[3*e + 0];
        double y  = rec_coords_[3*e + 1];
        double z  = rec_coords_[3*e + 2];

        // traction components
        double tx = traction[bs*e + 0];
        double ty = traction[bs*e + 1];
        double tz = traction[bs*e + 2];

        // write one line, then newline
        out
          << x << ',' << y << ',' << z << ','
          << tx << ',' << ty << ',' << tz << '\n';
    }

    out.close();
}
arma::vec multiply(const arma::vec &x) {
    // perm1_ = source permutation, perm0_ = receiver permutation
    std::size_t n_src = perm1_.size();
    std::size_t n_rec = perm0_.size();
    std::size_t n_in  = n_src * bs_;
    std::size_t n_out = n_rec * bs_;
    assert(x.n_elem == static_cast<long>(n_in));
    // 1) Permute input x into cluster order xc
    std::vector<double> xc(n_in, 0.0);
    std::cout << "xc.size() = " << xc.size() << std::endl;
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < n_src; ++i) {
        for (int d = 0; d < bs_; ++d) {
            xc[i * bs_ + d] = x[perm1_[i] * bs_ + d];  // source permutation
        }
    }
    // 2) Prepare output scratch yc
    std::vector<double> yc(n_out, 0.0);

    // 3) Apply full-rank blocks
    std::cout << "Applying full-rank blocks..." << std::endl;
    // #pragma omp parallel for schedule(dynamic)
    std::cout << "Number of full-rank blocks: " << fullRankBlocks_.size() << std::endl;
    for (size_t b = 0; b < fullRankBlocks_.size(); ++b) {
    const auto &blk = frBlocks_[b];
    const auto &A   = fullRankBlocks_[b];
    int i0 = blk.i0, i1 = blk.i1;
    int j0 = blk.j0, j1 = blk.j1;
    int ni = bs_ * (i1 - i0);
    int nj = bs_ * (j1 - j0);

    // slice the correct piece of xc:
    arma::vec xblock(&xc[j0 * bs_], nj, /*copy=*/false, /*strict=*/true);

    // wrap your raw il::Array2D into an Armadillo matrix:
    arma::Mat<double> Ablock(const_cast<double*>(A->data()), ni, nj, false, true);

    // sanity check before multiply:
    assert(Ablock.n_cols == xblock.n_rows && 
           "full-rank block shape mismatch");

    // perform the mat-vec:
    arma::vec yblock = Ablock * xblock;  // ni×nj * nj×1 → ni×1

    // scatter back into yc starting at receiver-cluster i0:
    #pragma omp critical
    for (int ii = 0; ii < ni; ++ii) {
        yc[i0 * bs_ + ii] += yblock[ii];
    }
    }
    std::cout << "Full-rank blocks applied." << std::endl;
    // 4) Apply low-rank blocks
    #pragma omp parallel for schedule(dynamic)
    for (size_t b = 0; b < lowRankBlocks_.size(); ++b) {
    const auto &blk = lrBlocks_[b];
    const auto &L   = *lowRankBlocks_[b];

    int i0 = blk.i0, i1 = blk.i1;
    int j0 = blk.j0, j1 = blk.j1;
    int ni = bs_ * (i1 - i0);
    int nj = bs_ * (j1 - j0);
    int r  = L.A.size(1);

    // 1) slice input
    arma::vec xblock(&xc[j0 * bs_], nj, /*copy=*/false, /*strict=*/true);

    // 2a) Bᵣⱼ * xblock  → tmp (r×1)
    arma::Mat<double> Bmat(const_cast<double*>(L.B.data()), r, nj, false, true);
    assert(Bmat.n_cols == xblock.n_rows);
    arma::vec tmp = Bmat * xblock;

    // 2b) Aᵢᵣ * tmp → yblock (ni×1)
    arma::Mat<double> Amat(const_cast<double*>(L.A.data()), ni, r, false, true);
    assert(Amat.n_cols == tmp.n_rows);
    arma::vec yblock = Amat * tmp;

    // 3) scatter back
    #pragma omp critical
    for (int ii = 0; ii < ni; ++ii) {
      yc[i0 * bs_ + ii] += yblock[ii];
    }
    }
    std::cout << "Low-rank blocks applied." << std::endl;

    // 5) Permute output back to global receiver order
    arma::vec y(n_out, arma::fill::zeros);
    #pragma omp parallel for schedule(static)
    for (std::size_t j = 0; j < n_rec; ++j) {
        for (int d = 0; d < bs_; ++d) {
            y[perm0_[j] * bs_ + d] = yc[j * bs_ + d];  // receiver permutation
        }
    }

    return y;
}
};
ExportBackgroundStress::ExportBackgroundStress(const Mesh& mesh,double E,double nu)
  : impl_(std::make_unique<Impl>(mesh,E,nu)) {}

ExportBackgroundStress::~ExportBackgroundStress() = default;
void ExportBackgroundStress::ComputeAndExport(const std::vector<double> &slip, const std::string &filename) {
  impl_->ComputeAndExport(slip, filename);
}

arma::vec ExportBackgroundStress::multiply(const arma::vec& x) const {
  return impl_->multiply(x);
}

} // namespace EQSim
