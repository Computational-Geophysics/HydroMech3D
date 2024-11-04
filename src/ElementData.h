//
// Created by Federico Ciardo on 10.08.21. All rights reserved.
//

// Inclusion from InsideLibrary
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/norm.h>

#ifndef INC_3DEQSIM_SRC_ELEMENTDATA_H
#define INC_3DEQSIM_SRC_ELEMENTDATA_H

namespace EQSim {

class ElementData {
 private:
  double theta_;                // Angle of s1 with respect to e1 = (1,0,0)
  double a_;                    // Half-length of element on one direction
  double b_;                    // Half-length of element on the other direction
  double areaElt_;     // Area of the element
  il::Array<double> n_{3, 0.};  // Normal vector
  il::Array<double> s1_{3, 0.};  // Shear vector 1
  il::Array<double> s2_{3, 0.};  // Shear vector 2
  il::Array<double> centroid_{3, 0.};

 public:
  ElementData() = default;
  ElementData(il::Array2D<double> &CoorElt, il::Array<double> &CentroidElt) {
    a_ = (std::sqrt(pow(fabs(CoorElt(1, 0) - CoorElt(2, 0)), 2) +
                    pow(fabs(CoorElt(1, 1) - CoorElt(2, 1)), 2) +
                    pow(fabs(CoorElt(1, 2) - CoorElt(2, 2)), 2))) /
         2.;

    b_ = (std::sqrt(pow(fabs(CoorElt(0, 0) - CoorElt(1, 0)), 2) +
                    pow(fabs(CoorElt(0, 1) - CoorElt(1, 1)), 2) +
                    pow(fabs(CoorElt(0, 2) - CoorElt(1, 2)), 2))) /
         2.;

    // The sign of the othonormal vectors follow the Right Hand rule
    il::Array2D<double> DiffCoorElt{3, 3, 0.};
    for (il::int_t I = 0; I < DiffCoorElt.size(0); ++I) {
      DiffCoorElt(I, 0) = -1. * CoorElt(I, 0) + CoorElt(I + 1, 0);
      DiffCoorElt(I, 1) = -1. * CoorElt(I, 1) + CoorElt(I + 1, 1);
      DiffCoorElt(I, 2) = -1. * CoorElt(I, 2) + CoorElt(I + 1, 2);
    }

    il::Array<double> aux{3, 0.};
    double normaux;
    aux[0] = (-1. * DiffCoorElt(0, 2) * DiffCoorElt(1, 1)) +
             DiffCoorElt(0, 1) * DiffCoorElt(1, 2);
    aux[1] = DiffCoorElt(0, 2) * DiffCoorElt(1, 0) -
             (DiffCoorElt(0, 0) * DiffCoorElt(1, 2));
    aux[2] = (-1. * DiffCoorElt(0, 1) * DiffCoorElt(1, 0)) +
             DiffCoorElt(0, 0) * DiffCoorElt(1, 1);
    normaux = il::norm(aux, il::Norm::L1);

    n_[0] = aux[0] / normaux;
    n_[1] = aux[1] / normaux;
    n_[2] = aux[2] / normaux;

    il::Array<double> aux1{3, 0.};
    double normaux1;
    aux1[0] = (-1. * (-1. * n_[2]) * DiffCoorElt(0, 1)) +
              (-1. * n_[1]) * DiffCoorElt(0, 2);
    aux1[1] =
        (-1. * n_[2]) * DiffCoorElt(0, 0) - ((-1. * n_[0]) * DiffCoorElt(0, 2));
    aux1[2] = (-1. * (-1. * n_[1]) * DiffCoorElt(0, 0)) +
              (-1. * n_[0]) * DiffCoorElt(0, 1);
    normaux1 = il::norm(aux1, il::Norm::L1);

    s1_[0] = aux1[0] / normaux1;
    s1_[1] = aux1[1] / normaux1;
    s1_[2] = aux1[2] / normaux1;

    il::Array<double> aux2{3, 0.};
    double normaux2;
    aux2[0] = (-1. * (-1. * n_[2]) * DiffCoorElt(1, 1)) +
              (-1. * n_[1]) * DiffCoorElt(1, 2);
    aux2[1] =
        (-1. * n_[2]) * DiffCoorElt(1, 0) - ((-1. * n_[0]) * DiffCoorElt(1, 2));
    aux2[2] = (-1. * (-1. * n_[1]) * DiffCoorElt(1, 0)) +
              (-1. * n_[0]) * DiffCoorElt(1, 1);
    normaux2 = il::norm(aux2, il::Norm::L1);

    s2_[0] = aux2[0] / normaux2;
    s2_[1] = aux2[1] / normaux2;
    s2_[2] = aux2[2] / normaux2;

    il::Array<double> e1{3, 0.};
    e1[0] = 1.;
    e1[1] = 0;
    e1[2] = 0;
    double costheta = (e1[0] * s1_[0] + e1[1] * s1_[1] + e1[2] * s1_[2]) /
               (il::norm(e1, il::Norm::L1) * il::norm(s1_, il::Norm::L1));
    theta_ = std::acos(costheta);
    if (s1_[1]<0){ theta_ = -1.* theta_; }

    centroid_ = CentroidElt;
    areaElt_ = ((2*a_) * (2*b_));
  }

  // Getter methods
  il::Array<double> getN() { return n_; };
  il::Array<double> getS1() { return s1_; };
  il::Array<double> getS2() { return s2_; };
  double getTheta() const { return theta_; };
  double getAreaElt() const { return areaElt_; };
  il::Array<double> getCentroidElt() const { return centroid_; };
  double getCentroidElt(il::int_t i) const { return centroid_[i]; };
  double getA() const { return a_; };
  double getB() const { return b_; };

};

}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_ELEMENTDATA_H
