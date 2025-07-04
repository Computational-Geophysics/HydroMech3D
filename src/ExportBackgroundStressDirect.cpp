#include "ExportBackgroundStressDirect.h"
#include "Mesh.h"

#include <armadillo>
#include <fstream>
#include <memory>
#include <array>
#include <iostream>

#include "core/be_mesh.h"
#include "elements/rectangle.h"
#include "elasticity/bie_elastostatic.h"
#include "hmat/hierarchical_representation.h"
#include "hmat/square_matrix_generator.h"
#include "io/bigwham_io_helper.h"
#include "hmat/bie_matrix_generator.h"

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
        std::cerr << "[ExportBackgroundStressDirect] Warning: "
                  << axis << "-panel normal not axis aligned\n";
    return sign * axisVal < 0;
}
}

namespace EQSim {

struct ExportBackgroundStressDirect::Impl {
    std::vector<double> rec_coords_;
    arma::Mat<double> dense_mat_;
    int bs_ = 3;

    std::shared_ptr<bigwham::Mesh> mesh_src_;
    std::shared_ptr<bigwham::Mesh> mesh_rec_;
    std::shared_ptr<bigwham::BieKernel<double>> ker_obj_;

    Impl(const Mesh& faultMesh, double E, double nu) {
        arma::mat coordsMat = faultMesh.getCoordinates();
        arma::imat connMat  = faultMesh.getNodesConnettivity();
        int nodesPerElt = connMat.n_cols;

        using src_elem = bigwham::Rectangle<0>;
        using rec_elem = bigwham::Rectangle<0>;
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
            for (int n = 0; n < nodesPerElt; ++n)
                conn_src.push_back(connMat(e,n));
        }

        // simple receiver mesh around origin similar to ExportBackgroundStress
        /*
        std::vector<double> zPlanes = { 0.0 };           // sample at the fault plane z=0
        int NX_ = 81, NY_ = 81;
        double x_min=-18., x_max=18., y_min=-18., y_max=18.;
        double dx = (x_max-x_min)/(NX_-1), dy=(y_max-y_min)/(NY_-1);
        double eps = 0.25 * std::min(dx,dy);

        std::vector<double> coor_rec;
        std::vector<int>    conn_rec;
        std::vector<double> centroids;
        int nodeIndexOffset = 0;

        for (double Z : zPlanes) {
          for (int i=0; i<NX_; ++i) {
            double X = x_min + i*dx;
            for (int j=0; j<NY_; ++j) {
              double Y = y_min + j*dy;
              // record the centroid
              centroids.push_back(X);
              centroids.push_back(Y);
              centroids.push_back(Z);
            
              // build a small quad in the XY-plane
              std::array<double,12> verts = {
                X+eps, Y+eps, Z,
                X-eps, Y+eps, Z,
                X-eps, Y-eps, Z,
                X+eps, Y-eps, Z
              };
          
              // compute its normal (should be ±Z)
              Vec3 n = normal(&verts[0], &verts[3], &verts[6]);
              int signWanted = +1;   // enforce normal pointing +Z
              if (needsFlip(n, signWanted, 'z')) {
                // reverse order so normal is +Z
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
              coor_rec.insert(coor_rec.end(), verts.begin(), verts.end());
            }
          }
        }
        rec_coords_ = centroids;*/
        std::vector<double> xPlanes = {0.0};
        int NY_=81;
        int NZ_=81;
        double y_min=-18.0, y_max=18.0;
        double z_min=-18.0, z_max=18.0;
        double dy=(y_max-y_min)/(NY_-1);
        double dz=(z_max-z_min)/(NZ_-1);
        double eps=0.25*std::min(dy,dz);
        std::vector<double> coor_rec;
        std::vector<int> conn_rec;
        std::vector<double> centroids;
        int nodeIndexOffset=0;
        for(double X:xPlanes){
            for(int j=0;j<NY_;++j){
                double Y=y_min+j*dy;
                for(int k=0;k<NZ_;++k){
                    double Z=z_min+k*dz;
                    centroids.push_back(X); centroids.push_back(Y); centroids.push_back(Z);
                    std::array<double,12> verts={
                        X, Y+eps, Z+eps,
                        X, Y-eps, Z+eps,
                        X, Y-eps, Z-eps,
                        X, Y+eps, Z-eps
                    };
                    Vec3 n = normal(&verts[0], &verts[3], &verts[6]);
                    coor_rec.insert(coor_rec.end(), verts.begin(), verts.end());
                    int signWanted = (X < 0.0 ? -1 : 1);
                    if(needsFlip(n, signWanted, 'x')){
                        conn_rec.insert(conn_rec.end(), {nodeIndexOffset+0,nodeIndexOffset+3,nodeIndexOffset+2,nodeIndexOffset+1});
                    }else{
                        conn_rec.insert(conn_rec.end(), {nodeIndexOffset+0,nodeIndexOffset+1,nodeIndexOffset+2,nodeIndexOffset+3});
                    }
                    nodeIndexOffset+=4;
                }
            }
        }
        rec_coords_ = centroids;
        bigwham::ElasticProperties elas(E,nu);
        mesh_src_ = bigwham::CreateMeshFromVect<src_elem>(3,4,coor_src,conn_src);
        mesh_rec_ = bigwham::CreateMeshFromVect<rec_elem>(3,4,coor_rec,conn_rec);
        ker_obj_ = std::make_shared<bigwham::BieElastostatic<src_elem,rec_elem,bigwham::ElasticKernelType::H>>(elas,3);

        // representation just for access to generator
        auto hr = HRepresentationRectangularMatrix(mesh_src_, mesh_rec_, 6561*3, 1e9);
        bigwham::BieMatrixGenerator<double> M(mesh_src_, mesh_rec_, ker_obj_, hr);

        il::Array2D<double> full(bs_ * mesh_rec_->num_elements(), bs_ * mesh_src_->num_elements());
        M.set(0,0,il::io, full.Edit());

        dense_mat_.set_size(full.size(0), full.size(1));
        for(il::int_t i=0;i<full.size(0);++i)
            for(il::int_t j=0;j<full.size(1);++j)
                dense_mat_(i,j) = full(i,j);
    }

    arma::vec multiply(const arma::vec& x) const {
        return dense_mat_ * x;
    }

    void ComputeAndExport(const std::vector<double>& slip, const std::string& filename){
        arma::vec slip_vec(slip.size());
        for(size_t i=0;i<slip.size();++i) slip_vec(i)=slip[i];
        arma::vec traction = multiply(slip_vec);
        std::ofstream out(filename);
        out << "x,y,z,tx,ty,tz\n";
        int N = rec_coords_.size()/3;
        for(int e=0;e<N;++e){
            out << rec_coords_[3*e] << ',' << rec_coords_[3*e+1] << ',' << rec_coords_[3*e+2] << ','
                << traction[bs_*e] << ',' << traction[bs_*e+1] << ',' << traction[bs_*e+2] << '\n';
        }
    }
};

ExportBackgroundStressDirect::ExportBackgroundStressDirect(const Mesh& mesh,double E,double nu)
    : impl_(std::make_unique<Impl>(mesh,E,nu)) {}

ExportBackgroundStressDirect::~ExportBackgroundStressDirect() = default;

void ExportBackgroundStressDirect::ComputeAndExport(const std::vector<double>& slip,const std::string& filename){
    impl_->ComputeAndExport(slip, filename);
}

} // namespace EQSim

