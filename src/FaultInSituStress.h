//
// Created by Federico Ciardo on 10.08.21. All rights reserved.
//

// Import from the project
#include "Mesh.h"

#ifndef INC_3DEQSIM_SRC_INSITUSTRESS_H
#define INC_3DEQSIM_SRC_INSITUSTRESS_H

namespace EQSim {

class FaultInSituStress {
 private:
  il::Array<double> sxx_;  // sigma_xx
  il::Array<double> syy_;  // sigma_yy
  il::Array<double> szz_;
  il::Array<double> sxy_;
  il::Array<double> sxz_;
  il::Array<double> syz_;

 public:
  // Constructors
  FaultInSituStress() = default;

  FaultInSituStress(il::Array<double> &sigmaxx, il::Array<double> &sigma_yy,
                    il::Array<double> &sigma_zz, il::Array<double> &sigma_xy,
                    il::Array<double> &sigma_xz, il::Array<double> &sigma_yz) {
    sxx_ = sigmaxx;
    syy_ = sigma_yy;
    szz_ = sigma_zz;
    sxy_ = sigma_xy;
    sxz_ = sigma_xz;
    syz_ = sigma_yz;
  };

  // Getter methods
  il::Array<double> getSxx() const { return sxx_; };
  il::Array<double> getSyy() const { return syy_; };
  il::Array<double> getSzz() const { return szz_; };
  il::Array<double> getSxy() const { return sxy_; };
  il::Array<double> getSxz() const { return sxz_; };
  il::Array<double> getSyz() const { return syz_; };
  double getSxx(il::int_t &i) const { return sxx_[i]; };
  double getSyy(il::int_t &i) const { return syy_[i]; };
  double getSzz(il::int_t &i) const { return szz_[i]; };
  double getSxy(il::int_t &i) const { return sxy_[i]; };
  double getSxz(il::int_t &i) const { return sxz_[i]; };
  double getSyz(il::int_t &i) const { return syz_[i]; };

  // Methods
  il::Array<double> InsituTractionsElt(EQSim::ElementData &my_elt,
                                       il::int_t &elt_idx) const {
    il::Array<double> n = my_elt.getN();
    il::Array<double> s1 = my_elt.getS1();
    il::Array<double> s2 = my_elt.getS2();

    il::Array<double> traction{3, 0.};

    // s1.Sig.n
    traction[0] = (s1[0] * sxx_[elt_idx] * n[0] + s1[0] * sxy_[elt_idx] * n[1] +
                   s1[0] * sxz_[elt_idx] * n[2]) +
                  (s1[1] * sxy_[elt_idx] * n[0] + s1[1] * syy_[elt_idx] * n[1] +
                   s1[1] * syz_[elt_idx] * n[2]) +
                  (s1[2] * sxz_[elt_idx] * n[0] + s1[2] * syz_[elt_idx] * n[1] +
                   s1[2] * szz_[elt_idx] * n[2]);  // ts1

    // s2.Sig.n
    traction[1] = (s2[0] * sxx_[elt_idx] * n[0] + s2[0] * sxy_[elt_idx] * n[1] +
                   s2[0] * sxz_[elt_idx] * n[2]) +
                  (s2[1] * sxy_[elt_idx] * n[0] + s2[1] * syy_[elt_idx] * n[1] +
                   s2[1] * syz_[elt_idx] * n[2]) +
                  (s2[2] * sxz_[elt_idx] * n[0] + s2[2] * syz_[elt_idx] * n[1] +
                   s2[2] * szz_[elt_idx] * n[2]);  // ts2

    // n.Sig.n
    traction[2] = (n[0] * sxx_[elt_idx] * n[0] + n[0] * sxy_[elt_idx] * n[1] +
                   n[0] * sxz_[elt_idx] * n[2]) +
                  (n[1] * sxy_[elt_idx] * n[0] + n[1] * syy_[elt_idx] * n[1] +
                   n[1] * syz_[elt_idx] * n[2]) +
                  (n[2] * sxz_[elt_idx] * n[0] + n[2] * syz_[elt_idx] * n[1] +
                   n[2] * szz_[elt_idx] * n[2]);  // tn

    return traction;
  };

  il::Array<double> AllInSituTractions(EQSim::Mesh &my_mesh) const {
    il::Array<double> inSituTractionsEltE{};
    il::Array<double> all_tractions{my_mesh.getNumberOfDofs(), 0.};

    // Number of Dofs per element
    il::int_t NDofsPerElt = 3;
    // Number of total elements
    il::int_t Nelts = my_mesh.getNumberOfElts();

    // Loop over all the elements
    for (il::int_t e = 0; e < Nelts; e++) {
      // Get data of element e
      EQSim::ElementData my_elt = my_mesh.getElementData(e);
      // Get tractions on element e
      inSituTractionsEltE = InsituTractionsElt(my_elt, e);
      // Fill the vector
      all_tractions[e * NDofsPerElt] = inSituTractionsEltE[0];
      all_tractions[e * NDofsPerElt + 1] = inSituTractionsEltE[1];
      all_tractions[e * NDofsPerElt + 2] = inSituTractionsEltE[2];
    }

    return all_tractions;
  }
};
}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_INSITUSTRESS_H
