//
// Created by Federico Ciardo on 01.10.21. All rights reserved.
//

// Inclusion from Standard library
#include <iostream>

// Import from the project
#include <JSON/nlohmann/json.hpp>

#include "ElementData.h"
#include "FluidProperties.h"
#include "Mesh.h"
#include "il/Array.h"

#ifndef INC_3DEQSIM_SRC_INJECTION_H
#define INC_3DEQSIM_SRC_INJECTION_H

// Using declaration for json library to lighten the code
using json = nlohmann::json;

namespace EQSim {

// This class encapsulates all the injection parameters
// Class members:
// - mass rates [Kg/s]
// - constant overpressure [Pa]
// - time at which injection change [s]
// - source elements, i.e. elements in which fluid is injected into
// - type of injection

class Injection {
 private:
  il::Array<double> mass_rates_source_elm_;
  il::Array<double> constant_overpressures_source_elm_;
  il::Array<double> time_at_which_injection_change_source_elm_;
  il::int_t source_element_;

  std::string type_of_injection_;

 public:
  // Constructor
  Injection() = default;

  Injection(json &j_injection) {
    if (j_injection.count("Type") == 1) {
      type_of_injection_ = j_injection["Type"].get<std::string>();
    } else {
      std::cout << "No injection type in input data!";
      il::abort();
    }

    if (j_injection.count("Source element") == 1) {
      source_element_ = j_injection["Source element"].get<il::int_t>();
    } else {
      std::cout << "No Source element in input data!";
      il::abort();
    }

    if (j_injection.count("Time at which injection change - Source Elem") ==
        1) {
      time_at_which_injection_change_source_elm_.Resize(
          j_injection["Time at which injection change - Source Elem"].size());
      for (il::int_t I = 0;
           I < time_at_which_injection_change_source_elm_.size(); ++I) {
        time_at_which_injection_change_source_elm_[I] =
            j_injection["Time at which injection change - Source Elem"][I];
      }
    } else {
      std::cout << "No Time at which injection change - Source Elem in "
                   "input data!";
      il::abort();
    }

    if (type_of_injection_ == "Rate control") {
      if (j_injection.count("Mass rate - Source Element") == 1) {
        mass_rates_source_elm_.Resize(
            j_injection["Mass rate - Source Element"].size());
        for (il::int_t I = 0; I < mass_rates_source_elm_.size(); ++I) {
          mass_rates_source_elm_[I] =
              j_injection["Mass rate - Source Element"][I];
        }
      } else {
        std::cout << "No Mass rate - Source Element in input data!";
        il::abort();
      }

      IL_ASSERT(time_at_which_injection_change_source_elm_.size() ==
                mass_rates_source_elm_.size());

    } else if (type_of_injection_ == "Pressure control") {
      if (j_injection.count("Constant overpressure - Source Element") == 1) {
        constant_overpressures_source_elm_.Resize(
            j_injection["Constant overpressure - Source Element"].size());
        for (il::int_t I = 0; I < constant_overpressures_source_elm_.size();
             ++I) {
          constant_overpressures_source_elm_[I] =
              j_injection["Constant overpressure - Source Element"][I];
        }
      } else {
        std::cout
            << "No constant overpressure - Source Element in input data!";
        il::abort();
      }

      IL_ASSERT(time_at_which_injection_change_source_elm_.size() ==
                constant_overpressures_source_elm_.size());
    }
  }

  // Getter methods
  il::Array<double> getTimeAtWhichInjectionChangeSourceElem() {
    return time_at_which_injection_change_source_elm_;
  };
  double getTimeAtWhichInjectionChangeSourceElem(il::int_t &i) {
    return time_at_which_injection_change_source_elm_[i];
  };
  il::int_t getSourceElement1() const { return source_element_; };
  il::Array<double> getMassRatesSourceElem() {
    return mass_rates_source_elm_;
  };
  il::Array<double> getConstantOverpressuresSourceElem() {
    return constant_overpressures_source_elm_;
  };

  // Setter methods
  void setTimeAtWhichInjectionChangeSourceElem(il::int_t &i, double &time) {
    time_at_which_injection_change_source_elm_[i] = time;
  };
  void setMassRateSourceElem(il::int_t &i, double &mass_rate) {
    mass_rates_source_elm_[i] = mass_rate;
  };
  void setConstantOverpressureSourceElem(il::int_t &i, double &const_dp) {
    constant_overpressures_source_elm_[i] = const_dp;
  };

  // Methods
  il::Array<double> calculate_Q_vector(EQSim::Mesh &Mesh,
                                       EQSim::FluidProperties &FluidProperties,
                                       double &curr_time) {
    il::int_t Nelts = Mesh.getNumberOfElts();
    il::Array<double> Q_vector{Nelts, 0.};
    il::int_t source_elt;
    il::Array<double> new_time_at_which_injection_change{};
    il::Array<double> new_mass_rate{};

    // If the type of injection is rate control, then calculates the Q vector
    // with the proper mass rates, otherwise keep the Q vector null
    if (type_of_injection_ == "Rate control") {
      if ((time_at_which_injection_change_source_elm_.size() >= 2) &&
          (curr_time >= time_at_which_injection_change_source_elm_[0])) {
        new_time_at_which_injection_change.Resize(
            time_at_which_injection_change_source_elm_.size() - 1);
        new_mass_rate.Resize(
            time_at_which_injection_change_source_elm_.size() - 1);
        for (il::int_t I = 0; I < new_time_at_which_injection_change.size();
             ++I) {
          new_time_at_which_injection_change[I] =
              time_at_which_injection_change_source_elm_[I + 1];
          new_mass_rate[I] = mass_rates_source_elm_[I + 1];
        }

        time_at_which_injection_change_source_elm_ =
            new_time_at_which_injection_change;
        mass_rates_source_elm_ = new_mass_rate;
      }

      EQSim::ElementData DataSourceElem1 =
          Mesh.getElementData(source_element_);
      for (il::int_t I = 0; I < mass_rates_source_elm_.size(); ++I) {
        Q_vector[source_element_] = (mass_rates_source_elm_[I] /
                                      (FluidProperties.getFluidDensity() *
                                       DataSourceElem1.getAreaElt()));  // m / s
      }
    }

    return Q_vector;
  }
};

}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_INJECTION_H
