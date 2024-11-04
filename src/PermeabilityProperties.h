//
// Created by Federico Ciardo on 30.08.21. All rights reserved.
//

// Inclusion from Standard library
#include <iostream>

// Import from the project
#include <JSON/nlohmann/json.hpp>

#include "Mesh.h"

#ifndef INC_3DEQSIM_SRC_PERMEABILITY_CUH
#define INC_3DEQSIM_SRC_PERMEABILITY_CUH

// Using declaration for json library to lighten the code
using json = nlohmann::json;

namespace EQSim {

// This class encapsulates all the fault permeability properties
// It includes the virtual function (interface) that sets the permeability
// parameters of a given perm law to each element of the mesh (based on the
// local MatID)

class PermeabilityProperties {
 protected:
  il::Array<double> permeabilities_;

 public:
  // Constructor
  PermeabilityProperties() = default;

  // Interfaces
  virtual void setPermeabilityParameters(json &j_perm_param,
                                         EQSim::Mesh &Mesh) = 0;

  // Getter functions
  il::Array<double> getPermeabilities() const { return permeabilities_; };
  double getPermeabilities(il::int_t i) const { return permeabilities_[i]; };
};

// This sub-class inherits from PermeabilityProperties class and encapsulates
// all the parameters of the constant permeability law. Since it does not have
// input parameters, we can get immediately the value of the permeability based
// on local matID

class ConstantPermeability : public PermeabilityProperties {
 public:
  void setPermeabilityParameters(json &j_perm_param,
                                 EQSim::Mesh &Mesh) override {
    // Initial check of keywords
    if (j_perm_param.count("Permeability") != 1) {
      std::cout << "Permeability keyword is wrong in input file! " << std::endl;
      il::abort();
    }

    il::Array<double> perm{1, 0.};
    IL_ASSERT(perm.size() == j_perm_param["Permeability"].size());

    for (il::int_t i = 0; i < 1; ++i) {
      perm[i] = j_perm_param["Permeability"][i];
    }

    // Assign all the material parameters to each element in the mesh

    il::int_t Nelts = Mesh.getNumberOfElts();
    il::Array<double> permeabilities{Nelts, 0.};
    il::int_t matID_idx;

    for (il::int_t I = 0; I < Nelts; ++I) {
      permeabilities[I] = perm[0];
    }

    permeabilities_ = permeabilities;
  };
};

}  // namespace EQSim
#endif  // INC_3DEQSIM_SRC_PERMEABILITY_CUH
