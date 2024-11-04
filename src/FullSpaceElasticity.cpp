//
// Created by Federico Ciardo on 11.08.21. All rights reserved.
//

// Inclusion from the project
#include "FullSpaceElasticity.h"

#include "StressKernelsP0/StressKernelsDxP0.h"
#include "StressKernelsP0/StressKernelsDyP0.h"
#include "StressKernelsP0/StressKernelsDzP0.h"
#include "Utils.cpp"

namespace EQSim {

il::Array2D<double> TractionsDueToDDsOnSingleEltP0(
    double &a, double &b, il::Array<double> &coor_centroid_elt,
    il::Array<double> &shear1_vector, il::Array<double> &shear2_vector,
    il::Array<double> &normal_vector, SolidMatrixProperties Matrix_Prop, il::io_t) {
  double nu = Matrix_Prop.getPoissonRatio();
  double Shear_Mod = Matrix_Prop.getShearModulus();

  double x1 = coor_centroid_elt[0];
  double x2 = coor_centroid_elt[1];
  double x3 = coor_centroid_elt[2];

  il::Array2D<double> StressTensorDueToDx{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDy{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDz{3, 3, 0.};

  // Stress tensor at (x1,x2,x3) due to Dx -> calculated via elastic kernel call
  StressTensorDueToDx =
      EQSim::StressTensorDueToDDxOnSingleEltP0(a, b, x1, x2, x3, nu, Shear_Mod);
  // Stress tensor at (x1,x2,x3) due to Dy -> calculated via elastic kernel call
  StressTensorDueToDy =
      EQSim::StressTensorDueToDDyOnSingleEltP0(a, b, x1, x2, x3, nu, Shear_Mod);
  // Stress tensor at (x1,x2,x3) due to Dz -> calculated via elastic kernel call
  StressTensorDueToDz =
      EQSim::StressTensorDueToDDzOnSingleEltP0(a, b, x1, x2, x3, nu, Shear_Mod);

  // Traction vector on (x1,x2,x3) due to Dx
  il::Array<double> TractionVectorDx =
      il::dot(StressTensorDueToDx, normal_vector);
  // Traction vector on (x1,x2,x3) due to Dy
  il::Array<double> TractionVectorDy =
      il::dot(StressTensorDueToDy, normal_vector);
  // Traction vector on (x1,x2,x3) due to Dz
  il::Array<double> TractionVectorDz =
      il::dot(StressTensorDueToDz, normal_vector);

  double ts1Dx = 0., ts1Dy = 0., ts1Dz = 0.;
  double ts2Dx = 0., ts2Dy = 0., ts2Dz = 0.;
  double tnDx = 0., tnDy = 0., tnDz = 0.;

  for (int I = 0; I < TractionVectorDx.size(); ++I) {
    ts1Dx += shear1_vector[I] * TractionVectorDx[I];
    ts1Dy += shear1_vector[I] * TractionVectorDy[I];
    ts1Dz += shear1_vector[I] * TractionVectorDz[I];

    ts2Dx += shear2_vector[I] * TractionVectorDx[I];
    ts2Dy += shear2_vector[I] * TractionVectorDy[I];
    ts2Dy += shear2_vector[I] * TractionVectorDz[I];

    tnDx += normal_vector[I] * TractionVectorDx[I];
    tnDy += normal_vector[I] * TractionVectorDy[I];
    tnDz += normal_vector[I] * TractionVectorDz[I];
  }

  il::Array2D<double> TractionsDueToDDsOnSingleEltP0{3, 3, 0.};

  TractionsDueToDDsOnSingleEltP0(0, 0) = ts1Dx;
  TractionsDueToDDsOnSingleEltP0(0, 1) = ts1Dy;
  TractionsDueToDDsOnSingleEltP0(0, 2) = ts1Dz;

  TractionsDueToDDsOnSingleEltP0(1, 0) = ts2Dx;
  TractionsDueToDDsOnSingleEltP0(1, 1) = ts2Dy;
  TractionsDueToDDsOnSingleEltP0(1, 2) = ts2Dz;

  TractionsDueToDDsOnSingleEltP0(2, 0) = tnDx;
  TractionsDueToDDsOnSingleEltP0(2, 1) = tnDy;
  TractionsDueToDDsOnSingleEltP0(2, 2) = tnDz;

  return TractionsDueToDDsOnSingleEltP0;
}

}  // namespace EQSim