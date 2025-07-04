// Hmatrix_bigwham.cpp
#include "Hmatrix_bigwham.h"
#include "Mesh.h"                               // your real Mesh + ElementData

// BigWham internals
#include "core/be_mesh.h"
#include "elements/rectangle.h"
#include "core/elastic_properties.h"
#include "elasticity/fullspace_iso_3d_rectangle/bie_elastostatic_rectangle_0_influence.h"
#include "hmat/hierarchical_representation.h"
#include "hmat/square_matrix_generator.h"
#include "hmat/hmatrix/hmat.h"

#include <unordered_map>
#include <cstring>  // for memcpy
#include <omp.h>

namespace EQSim {

struct HMatrixBigwham::Impl {
  std::vector<std::unique_ptr<bigwham::LowRank<double>>> lowRankBlocks;
  std::vector<std::unique_ptr<il::Array2D<double>>>      fullRankBlocks;
  std::vector<int> perm0, perm1;
  struct FR { int i0,i1,j0,j1; };
  struct LR { int i0,i1,j0,j1; };
  std::vector<FR> frBlocks;
  std::vector<LR> lrBlocks;
  int bs;
  Impl(const Mesh& mesh, double E, double nu) {
    using Rec0   = bigwham::Rectangle<0>;
    using BMesh  = bigwham::BEMesh<Rec0>;
    using Kernel = bigwham::BieElastostatic<Rec0,Rec0,bigwham::ElasticKernelType::H>;
    using MatrixGen = bigwham::SquareMatrixGenerator<double>;
    omp_set_dynamic(1);        // disable dynamic adjustment
    omp_set_num_threads(8);    // hard‐code 8 threads

    arma::mat coordsMat = mesh.getCoordinates();          // Nnodes x 3
    arma::imat connMat   = mesh.getNodesConnettivity();  // Nelts   x 4

    std::vector<std::array<double,3>> coords;
    coords.reserve(coordsMat.n_rows);
    //#pragma omp parallel for schedule(static)
    for (arma::uword i = 0; i < coordsMat.n_rows; ++i) {
      coords.push_back({ coordsMat(i,0), coordsMat(i,1), coordsMat(i,2) });
    }
    std::vector<std::array<int,4>> conn;
    conn.reserve(connMat.n_rows);
    //#pragma omp parallel for schedule(static)
    for (arma::uword e = 0; e < connMat.n_rows; ++e) {
      conn.push_back({
        static_cast<int>(connMat(e,0)),
        static_cast<int>(connMat(e,1)),
        static_cast<int>(connMat(e,2)),
        static_cast<int>(connMat(e,3))
      });
    }

    // Build BigWham BEMesh 
    il::Array2D<double>    coordsArr(coords.size(),    3);
    il::Array2D<il::int_t> connArr(conn.size(),        4);
    for (std::size_t i = 0; i < coords.size(); ++i)
      for (int d = 0; d < 3; ++d)
        coordsArr(i,d) = coords[i][d];
    for (std::size_t e = 0; e < conn.size(); ++e)
      for (int k = 0; k < 4; ++k)
        connArr(e,k) = static_cast<il::int_t>(conn[e][k]);

    auto mesh_sp = std::make_shared<BMesh>(coordsArr, connArr);
    mesh_sp->ConstructMesh();

    // Kernel + clustering parameters
    bigwham::ElasticProperties elas(E, nu);
    auto kernel_sp = std::make_shared<Kernel>(elas, /*spaceDim=*/3);
    std::size_t leafSize = 64;
    double      eta      = 2.0;
    double      epsACA   = 1e-4;

    auto hr_sp = bigwham::HRepresentationSquareMatrix(mesh_sp, leafSize, eta);
    const auto& patt = hr_sp->pattern_;
    bs = 3; // always 3 DOFs per collocation

    // copy permutations
    {
      auto &p0 = hr_sp->permutation_0_, &p1 = hr_sp->permutation_1_;
      perm0.resize(p0.size());
      perm1.resize(p1.size());
      for (int i = 0; i < (int)p0.size(); ++i) {
        perm0[i] = p0[i];
        perm1[i] = p1[i];
      }
    }

    // Build blocks via ACA + full-rank extraction 
    MatrixGen gen(mesh_sp, kernel_sp, hr_sp);

    // low‐rank
    lowRankBlocks.reserve(patt.n_LRB);
    lrBlocks.reserve(patt.n_LRB);
    //#pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < patt.n_LRB; ++k) {
      int i0 = patt.LRB_pattern(1,k), i1 = patt.LRB_pattern(3,k);
      int j0 = patt.LRB_pattern(2,k), j1 = patt.LRB_pattern(4,k);
      auto lra = bigwham::adaptiveCrossApproximation<3>(
        gen,{i0,i1},{j0,j1},epsACA);
      lowRankBlocks.emplace_back(std::move(lra));
      lrBlocks.push_back({i0,i1,j0,j1});
    }

    // full‐rank
    fullRankBlocks.reserve(patt.n_FRB);
    frBlocks.reserve(patt.n_FRB);
    //#pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < patt.n_FRB; ++k) {
      int i0 = patt.FRB_pattern(1,k), i1 = patt.FRB_pattern(3,k);
      int j0 = patt.FRB_pattern(2,k), j1 = patt.FRB_pattern(4,k);
      int ni = bs * (i1 - i0), nj = bs * (j1 - j0);
      auto A = std::make_unique<il::Array2D<double>>(ni,nj);
      gen.set(i0,j0,il::io,A->Edit());
      fullRankBlocks.emplace_back(std::move(A));
      frBlocks.push_back({i0,i1,j0,j1});
    }

  }
  
  arma::vec multiply(const arma::vec& x) const {
    std::size_t Ncoll = perm1.size(), N = x.n_elem;
    std::vector<double> xc(N, 0.0), yc(N, 0.0);
  
    // 1) permute in
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < Ncoll; ++i)
      for (int j = 0; j < bs; ++j)
        xc[i*bs + j] = x[perm1[i]*bs + j];
  
    // 2) full-rank blocks via Armadillo
    #pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < (int)fullRankBlocks.size(); ++b) {
      auto& blk = frBlocks[b];
      auto& A   = *fullRankBlocks[b];
  
      int i0 = blk.i0, i1 = blk.i1;
      int j0 = blk.j0, j1 = blk.j1;
      int ni = bs*(i1 - i0), nj = bs*(j1 - j0);
  
      // view the relevant slice of xc
      arma::vec xblock(&xc[j0*bs], nj, false, true);
      arma::Mat<double> Ablock(const_cast<double*>(A.data()), ni, nj, false, true);

      // perform the mat-vec
      arma::vec yblock = Ablock * xblock;
  
      // scatter back
      #pragma omp critical
      for (int ii = 0; ii < ni; ++ii) {
        yc[(i0*bs) + ii] += yblock[ii];
      }
    }
  
    // low-rank blocks via Armadillo
    #pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < (int)lowRankBlocks.size(); ++b) {
      auto& blk = lrBlocks[b];
      auto& L   = *lowRankBlocks[b];
  
      int i0 = blk.i0, i1 = blk.i1;
      int j0 = blk.j0, j1 = blk.j1;
      int ni = bs*(i1 - i0), nj = bs*(j1 - j0);
      int r  = L.A.size(1);
  
      arma::vec xblock(&xc[j0*bs], nj, false, true);
  
      //arma::Mat<double> Bmat(r, nj);
      //std::memcpy(Bmat.memptr(), L.B.data(), sizeof(double)*r*nj);
      arma::Mat<double> Bmat(const_cast<double*>(L.B.data()), r, nj, false, true);
      arma::vec tmp = Bmat * xblock;
  
      //arma::Mat<double> Amat(ni, r);
      //std::memcpy(Amat.memptr(), L.A.data(), sizeof(double)*ni*r);
      arma::Mat<double> Amat(const_cast<double*>(L.A.data()), ni, r, false, true);
      arma::vec yblock = Amat * tmp;
  
      #pragma omp critical
      for (int ii = 0; ii < ni; ++ii) {
        yc[(i0*bs) + ii] += yblock[ii];
      }
    }
  
    // permute back out
    arma::vec y(N);
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < Ncoll; ++i) {
      for (int j = 0; j < bs; ++j) {
        y[perm0[i]*bs + j] = yc[i*bs + j];
      }
    }
  
    return y;
  }
  
};

HMatrixBigwham::HMatrixBigwham(const Mesh& mesh,double E,double nu)
  : impl_(std::make_unique<Impl>(mesh,E,nu)) {}

HMatrixBigwham::~HMatrixBigwham() = default;

arma::vec HMatrixBigwham::multiply(const arma::vec& x) const {
  return impl_->multiply(x);
}

} // namespace EQSim
