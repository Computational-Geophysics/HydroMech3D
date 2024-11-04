//
// Created by Federico Ciardo on 30.08.21. All rights reserved.
//

// Inclusion from Standard library
#include <iostream>

// Import from the project
#include "Mesh.h"
#include "Solution.h"

#ifndef INC_3DEQSIM_SRC_FRICTION_H
#define INC_3DEQSIM_SRC_FRICTION_H

namespace EQSim {

// This class encapsulates all the friction properties
// It includes the virtual function (interface) that sets the friction
// parameters of a given friction law to each element of the mesh (based on the
// local MatID)

class FrictionProperties {
 protected:
  il::Array<double> friction_coefficients_;
  il::Array<double> a_values_;
  il::Array<double> b_values_;
  il::Array<double> reference_velocities_;
  il::Array<double> state_evolution_distances_;
  il::Array<double> reference_friction_coefficients_;
  il::Array<double> peak_friction_coefficients_;
  il::Array<double> residual_friction_coefficients_;
  il::Array<double> residual_slips_;

 public:
  // Constructor
  FrictionProperties() = default;

  // Interface
  virtual void setFrictionParameters(json &j_friction_param,
                                     EQSim::Mesh &Mesh) = 0;

  // Getter methods
  il::Array<double> getFrictionCoefficient() { return friction_coefficients_; }
  double getFrictionCoefficient(il::int_t i) {
    return friction_coefficients_[i];
  }
  il::Array<double> get_a_values() const { return a_values_; };
  double get_a_values(il::int_t &i) const { return a_values_[i]; };
  il::Array<double> get_b_values() const { return b_values_; };
  double get_b_values(il::int_t &i) const { return b_values_[i]; };
  il::Array<double> get_reference_velocities() {
    return reference_velocities_;
  };
  double get_reference_velocities(il::int_t &i) {
    return reference_velocities_[i];
  };
  il::Array<double> get_state_evolution_distances_() const {
    return state_evolution_distances_;
  };
  double get_state_evolution_distances_(il::int_t &i) const {
    return state_evolution_distances_[i];
  };
  il::Array<double> get_reference_friction_coefficients_() const {
    return reference_friction_coefficients_;
  };
  double get_reference_friction_coefficients_(il::int_t i) const {
    return reference_friction_coefficients_[i];
  };
};

// This sub-class inherits from FrictionProperties class and encapsulates all
// the parameters of the Rate and State friction law
// The first five members are used to import the input parameters from json
// input file.
// The last five members, instead, include the input parameters for
// all the mesh.

class RateAndStateFriction : public FrictionProperties {
 private:
  il::Array<double> a_;
  il::Array<double> b_;
  il::Array<double> Vo_;
  il::Array<double> Dc_;
  il::Array<double> fo_;

 public:
  void setFrictionParameters(json &j_friction_param,
                             EQSim::Mesh &Mesh) override {
    // Initial check of keywords
    if (j_friction_param.count("Reference friction coefficient") != 1) {
      std::cout
          << "Reference friction coefficient keyword is wrong in input file! "
          << std::endl;
      il::abort();
    }

    if (j_friction_param.count("a-values") != 1) {
      std::cout << "a-values keyword is wrong in input file! " << std::endl;
      il::abort();
    }

    if (j_friction_param.count("b-values") != 1) {
      std::cout << "b-values keyword is wrong in input file! " << std::endl;
      il::abort();
    }

    if (j_friction_param.count("State evolution distances") != 1) {
      std::cout << "State evolution distances keyword is wrong in input file! "
                << std::endl;
      il::abort();
    }

    if (j_friction_param.count("Reference velocities") != 1) {
      std::cout << "Reference velocities keyword is wrong in input file! "
                << std::endl;
      il::abort();
    }

    il::Array<double> a{1, 0.};
    il::Array<double> b{1, 0.};
    il::Array<double> Dc{1, 0.};
    il::Array<double> Vo{1, 0.};
    il::Array<double> fo{1, 0.};

    IL_ASSERT(fo.size() ==
              j_friction_param["Reference friction coefficient"].size());
    IL_ASSERT(a.size() == j_friction_param["a-values"].size());
    IL_ASSERT(b.size() == j_friction_param["b-values"].size());
    IL_ASSERT(Dc.size() ==
              j_friction_param["State evolution distances"].size());
    IL_ASSERT(Vo.size() == j_friction_param["Reference velocities"].size());

    for (il::int_t i = 0; i < 1; ++i) {
      fo[i] = j_friction_param["Reference friction coefficient"][i];
      a[i] = j_friction_param["a-values"][i];
      b[i] = j_friction_param["b-values"][i];
      Dc[i] = j_friction_param["State evolution distances"][i];
      Vo[i] = j_friction_param["Reference velocities"][i];
    }

    fo_ = fo;
    a_ = a;
    b_ = b;
    Dc_ = Dc;
    Vo_ = Vo;

    // Assign all the material parameters to each element in the mesh

    il::int_t Nelts = Mesh.getNumberOfElts();
    il::Array<double> foValues{Nelts, 0.};
    il::Array<double> aValues{Nelts, 0.};
    il::Array<double> bValues{Nelts, 0.};
    il::Array<double> reference_velocities{Nelts, 0.};
    il::Array<double> state_evolution_distances{Nelts, 0.};

    for (il::int_t I = 0; I < Nelts; ++I) {
      foValues[I] = fo_[0];
      aValues[I] = a_[0];
      bValues[I] = b_[0];
      reference_velocities[I] = Vo_[0];
      state_evolution_distances[I] = Dc_[0];
    }

    reference_friction_coefficients_ = foValues;
    a_values_ = aValues;
    b_values_ = bValues;
    reference_velocities_ = reference_velocities;
    state_evolution_distances_ = state_evolution_distances;
  };
};

}  // namespace EQSim
#endif  // INC_3DEQSIM_SRC_FRICTION_H
