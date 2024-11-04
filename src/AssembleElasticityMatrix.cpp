//
// Created by Federico Ciardo on 11.10.22. All rights reserved.
//

// Import from the project
#include "AssembleElasticityMatrix.h"

namespace EQSim {

il::Array2D<double> RotationMatrix3D(il::Array<double> &normal_vector,
                                            double &theta) {
  il::Array2D<double> R{3, 3, 0.};

  double n1, n2, n3, n1square, n2square, n3square;
  n1 = normal_vector[0];
  n2 = normal_vector[1];
  n3 = normal_vector[2];
  n1square = n1 * n1;
  n2square = n2 * n2;
  n3square = n3 * n3;

  R(0, 0) = cos(theta) + (n1square * (1 - cos(theta)));
  R(0, 1) = (n1 * n2 * (1 - cos(theta))) - (n3 * sin(theta));
  R(0, 2) = (n1 * n3 * (1 - cos(theta))) + (n2 * sin(theta));

  R(1, 0) = (n1 * n2 * (1 - cos(theta))) + n3 * sin(theta);
  R(1, 1) = cos(theta) + (n2square * (1 - cos(theta)));
  R(1, 2) = (n2 * n3 * (1 - cos(theta))) - (n1 * sin(theta));

  R(2, 0) = (n1 * n3 * (1 - cos(theta))) - (n2 * sin(theta));
  R(2, 1) = (n2 * n3 * (1 - cos(theta))) + (n1 * sin(theta));
  R(2, 2) = cos(theta) + (n3square * (1 - cos(theta)));

  return R;
}

il::Array2D<double> AssembleElastMat(Mesh &mesh, EQSim::SolidMatrixProperties &matrix_prop) {

  il::Array<double> normal_source_elt{};
  il::Array<double> normal_receiver_elt{};
  il::Array<double> s1_receiver_elt{};
  il::Array<double> s2_receiver_elt{};
  il::Array<double> ne{};
  il::Array<double> se1{};
  il::Array<double> se2{};
  il::Array<double> relative_distance_receiv_source{3, 0.};
  il::Array<double> Xe{3, 0.};
  il::Array2D<double> TractionsDueToDDsOnSingleElt{3, 3, 0.};
  il::Array2D<double> R{3, 3, 0.};
  il::Array2D<double> Rt{3, 3, 0.};
  double theta_source_elt;
  EQSim::ElementData elt_s;
  EQSim::ElementData elt_r;

  il::Array2D<double> ElastMat{mesh.getNumberOfDofs(), mesh.getNumberOfDofs(), 0.};

  for (il::int_t e = 0; e < mesh.getNumberOfElts(); ++e) {
    elt_s = mesh.getElementData(e);
    for (il::int_t j = 0; j < mesh.getNumberOfElts(); ++j) {
      elt_r = mesh.getElementData(j);
      normal_source_elt = elt_s.getN();
      theta_source_elt = elt_s.getTheta();
      R = EQSim::RotationMatrix3D(normal_source_elt, theta_source_elt);
      Rt = R;
      Rt(0, 1) = R(1, 0);
      Rt(0, 2) = R(2, 0);
      Rt(1, 0) = R(0, 1);
      Rt(1, 2) = R(2, 1);
      Rt(2, 0) = R(0, 2);
      Rt(2, 1) = R(1, 2);

      for (il::int_t I = 0; I < relative_distance_receiv_source.size();
           ++I) {
        relative_distance_receiv_source[I] =
            elt_r.getCentroidElt(I) - elt_s.getCentroidElt(I);
      }
      Xe = il::dot(Rt, relative_distance_receiv_source);
      normal_receiver_elt = elt_r.getN();
      s1_receiver_elt = elt_r.getS1();
      s2_receiver_elt = elt_r.getS2();
      se1 = il::dot(Rt, s1_receiver_elt);
      se2 = il::dot(Rt, s2_receiver_elt);
      ne = il::dot(Rt, normal_receiver_elt);

      auto a_elt_s = elt_s.getA();
      auto b_elt_s = elt_s.getB();

      TractionsDueToDDsOnSingleElt =
          EQSim::TractionsDueToDDsOnSingleEltP0(
              a_elt_s, b_elt_s, Xe, se1, se2, ne, matrix_prop, il::io);

      for (il::int_t j2 = 0; j2 < TractionsDueToDDsOnSingleElt.size(0);
           ++j2) {
        for (il::int_t i2 = 0; i2 < TractionsDueToDDsOnSingleElt.size(1);
             ++i2) {
          ElastMat(j * 3 + i2, e * 3 + j2) = TractionsDueToDDsOnSingleElt(i2, j2);
        }
      }

    }
  }

  return ElastMat;
}


}