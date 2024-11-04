//
// Created by Federico Ciardo on 30.09.21. All rights reserved.
//

#include "il/Array.h"

#ifndef INC_3DEQSIM_SRC_FAULTPROPERTIES_H
#define INC_3DEQSIM_SRC_FAULTPROPERTIES_H

namespace EQSim {

class FrictionProperties;
class DilatancyProperties;
class PermeabilityProperties;

// This class encapsulates all the fault properties
// Class members:
// - fault gauge porosity [-]
// - void volume compressibility [1/Pa]
// - fault hydraulic aperture [m]

class FaultProperties {
 private:
  double porosity_;
  double void_compressibility_;
  il::Array<double> initial_hydraulic_aperture_;
  EQSim::PermeabilityProperties *permeabilityProperties_;
  EQSim::FrictionProperties *friction_properties_;
  il::Array<double> initial_DDs_;
  il::Array<double> initial_DDs_rates_;
  il::Array<double> initial_state_variables_;
  il::Array<double> initial_plastic_porosity_;
  il::Array<double> internal_loads_;

  // Constructors
 public:
  FaultProperties() = default;
  FaultProperties(double &porosity, double &void_compressibility,
                  il::Array<double> &initial_hydraulic_aperture,
                  EQSim::FrictionProperties *friction_properties,
                  EQSim::PermeabilityProperties *permeabilityProperties,
                  il::Array<double> &initial_DDs,
                  il::Array<double> &initial_DDs_rates,
                  il::Array<double> &initial_state_variables,
                  il::Array<double> &initial_plastic_porosity,
                  il::Array<double> &internal_loads) {
    porosity_ = porosity;
    void_compressibility_ = void_compressibility;
    initial_hydraulic_aperture_ = initial_hydraulic_aperture;
    friction_properties_ = friction_properties;
    permeabilityProperties_ = permeabilityProperties;
    initial_DDs_ = initial_DDs;
    initial_DDs_rates_ = initial_DDs_rates;
    initial_state_variables_ = initial_state_variables;
    initial_plastic_porosity_ = initial_plastic_porosity;
    internal_loads_ = internal_loads;
  }

  // Getter methods
  double getFaultPorosity() const { return porosity_; };
  double getInitialPlasticFaultPorosity(il::int_t index) const {
    return initial_plastic_porosity_[index];
  };
  double getFaultVoidCompressibility() const { return void_compressibility_; };
  double getInitialFaultHydraulicAperture(il::int_t i) const {
    return initial_hydraulic_aperture_[i];
  };
  PermeabilityProperties *getPermeabilityProperties() {
    return permeabilityProperties_;
  };
  FrictionProperties *getFrictionProperties() { return friction_properties_; };
  il::Array<double> getInitialDDs() { return initial_DDs_; };
  il::Array<double> getInitialDDsRates() { return initial_DDs_rates_; };
  il::Array<double> getInitialStateVariables() {
    return initial_state_variables_;
  };
  il::Array<double> getInternalLoads() { return internal_loads_; };
  il::Array<double> getAmbientPressureDistribution() {
    il::Array<double> ambient_pressure{internal_loads_.size() / 3, 0.};
    for (il::int_t I = 0; I < internal_loads_.size() / 3; ++I) {
      ambient_pressure[I] = internal_loads_[3 * I + 2];
    }
    return ambient_pressure;
  };

  // Setter methods
  void setFaultPorosity(double &porosity) { porosity_ = porosity; };
  void setFaultVoidCompressibility(double &void_compressibility) {
    void_compressibility_ = void_compressibility;
  };
};
}  // namespace EQSim
#endif  // INC_3DEQSIM_SRC_FAULTPROPERTIES_H
