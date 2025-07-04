#pragma once
#include <armadillo>
#include <memory>

namespace EQSim {
  class Mesh;
}

namespace EQSim {

class HMatrixBigwham {
public:
  /// Build the H-matrix for 'mesh' with Young’s modulus E and Poisson ratio nu.
  HMatrixBigwham(const EQSim::Mesh& mesh, double E, double nu);

  ~HMatrixBigwham();
  
  arma::vec multiply(const arma::vec& x) const;

  private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace EQSim

