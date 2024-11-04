//
// Created by Federico Ciardo on 30.08.21. All rights reserved.
//

// Import from IL library
#include <il/Array.h>
#include <il/Array2D.h>

#include <JSON/nlohmann/json.hpp>

#include "FaultInSituStress.h"
#include "FaultProperties.h"
#include "FluidProperties.h"
#include "Injection.h"
#include "Mesh.h"
#include "SolidMatrixProperties.h"
#include "SolverParameters.h"

// Using declaration for json library to lighten the code
using json = nlohmann::json;

#ifndef INC_3DEQSIM_SRC_SOLUTION_H
#define INC_3DEQSIM_SRC_SOLUTION_H

namespace EQSim {

// Forward declaration to avoid circular dependency between header files
class FrictionProperties;
class PermeabilityProperties;

class Solution {
 protected:
  double current_time_;
  double time_step_;
  EQSim::SolverParameters solver_parameters_;
  il::Array<double> insitu_tractions_;
  il::Array<double> total_tractions_;
  il::Array<double> DDs_;
  il::Array<double> DDs_rates_;
  il::Array<double> pressure_;
  il::Array<double> pressure_o_;
  EQSim::Injection injection_;
  il::Array<double> hydraulic_aperture_;
  EQSim::Mesh mesh_;

 public:
  Solution() = default;

  Solution(double &current_time, double &time_step,
           EQSim::SolverParameters &solver_parameters,
           il::Array<double> &insitu_tractions, il::Array<double> &DDs,
           il::Array<double> &DDs_rates, il::Array<double> &pressure,
           il::Array<double> &init_pressure,
           il::Array<double> &hydraulic_aperture, EQSim::Mesh &mesh) {
    mesh_ = mesh;
    current_time_ = current_time;
    time_step_ = time_step;
    solver_parameters_ = solver_parameters;
    insitu_tractions_ = insitu_tractions;
    DDs_rates_ = DDs_rates;
    DDs_ = DDs;
    pressure_ = pressure;
    pressure_o_ = init_pressure;
    hydraulic_aperture_ = hydraulic_aperture;
  }

  // Methods
  json createJsonObjectToExport();
  void exportBaseSolutionToJson(json &j_obj_to_export, std::string &filename);
  void exportBaseSolutionToUBJson(json &j_obj_to_export, std::string &filename);

  // Getter methods
  EQSim::Mesh getMesh() const { return mesh_; };
  EQSim::SolverParameters getSolverParametersStructure() {
    return solver_parameters_;
  };
  il::Array<double> getDDs() const { return DDs_; };
  double getDDs(il::int_t &i) const { return DDs_[i]; };
  il::Array<double> getSlip1() const {
    il::Array<double> slip1{DDs_.size() / 3, 0.};
    for (il::int_t I = 0; I < slip1.size(); ++I) {
      slip1[I] = DDs_[3 * I];
    }
    return slip1;
  };
  il::Array<double> getSlip2() const {
    il::Array<double> slip2{DDs_.size() / 3, 0.};
    for (il::int_t I = 0; I < slip2.size(); ++I) {
      slip2[I] = DDs_[3 * I + 1];
    }
    return slip2;
  };
  il::Array<double> getOpening() const {
    il::Array<double> opening{DDs_.size() / 3, 0.};
    for (il::int_t I = 0; I < opening.size(); ++I) {
      opening[I] = DDs_[3 * I + 2];
    }
    return opening;
  };
  il::Array<double> getDDsRates() const { return DDs_rates_; };
  double getDDsRates(il::int_t index) const { return DDs_rates_[index]; };
  il::Array<double> getSlipRate1() const {
    il::Array<double> slip_rate1{DDs_rates_.size() / 3, 0.};
    for (il::int_t I = 0; I < slip_rate1.size(); ++I) {
      slip_rate1[I] = DDs_rates_[3 * I];
    }
    return slip_rate1;
  };
  il::Array<double> getSlipRate2() const {
    il::Array<double> slip_rate2{DDs_rates_.size() / 3, 0.};
    for (il::int_t I = 0; I < slip_rate2.size(); ++I) {
      slip_rate2[I] = DDs_rates_[3 * I + 1];
    }
    return slip_rate2;
  };
  il::Array<double> getPressure() const { return pressure_; };
  double getPressure(il::int_t i) const { return pressure_[i]; };
  il::Array<double> getInitPressure() const { return pressure_o_; };
  double getInitPressure(il::int_t &i) const { return pressure_o_[i]; };
  il::Array<double> getInsituTractions() const { return insitu_tractions_; };
  double getInsituTractions(il::int_t &i) const {
    return insitu_tractions_[i];
  };
  il::Array<double> getNormalInSituTractions() const {
    il::Array<double> normal_tractions{insitu_tractions_.size() / 3, 0.};
    for (int I = 0; I < normal_tractions.size(); ++I) {
      normal_tractions[I] = insitu_tractions_[3 * I + 2];
    }
    return normal_tractions;
  };
  double getNormalInSituTractions(il::int_t &i) const {
    return insitu_tractions_[3 * i + 2];
  };
  double getNormalTotalTractions(il::int_t &i) const {
    return total_tractions_[3 * i + 2];
  };
  double getCurrentTime() const { return current_time_; };
  double getTimeStep() const { return time_step_; };
  il::Array<double> getHydraulicAperture() const {
    return hydraulic_aperture_;
  };
  double getHydraulicAperture(il::int_t &i) const {
    return hydraulic_aperture_[i];
  };
  double getSlipRateAlongSlipDirection(il::int_t &elm_i,
                                       il::int_t &slip_direction) {
    return DDs_rates_[3 * elm_i + slip_direction];
  }

  // Setter methods
  void setSlip1(double &slip1, il::int_t &i) { DDs_[3 * i] = slip1; }
  void setSlip2(double &slip2, il::int_t &i) { DDs_[3 * i + 1] = slip2; }
  void setOpening(double &opening, il::int_t &i) { DDs_[3 * i + 2] = opening; }
  void setSlipRate1(double &slip_rate1, il::int_t &i) {
    DDs_rates_[3 * i] = slip_rate1;
  }
  void setSlipRate2(double &slip_rate2, il::int_t &i) {
    DDs_rates_[3 * i + 1] = slip_rate2;
  }
  void setSlipAlongSlipDirection(il::Array<double> &slip,
                                 il::int_t &slip_direction) {
    for (il::int_t I = 0; I < DDs_.size() / 3; ++I) {
      DDs_[3 * I + slip_direction] = slip[I];
    }
  }
  void setSlipAlongSlipDirection(il::int_t &i, const double &slip,
                                 il::int_t &slip_direction) {
    DDs_[3 * i + slip_direction] = slip;
  }
  void setSlipRateAlongSlipDirection(il::Array<double> &slip_rate,
                                     il::int_t &slip_direction) {
    for (il::int_t I = 0; I < DDs_rates_.size() / 3; ++I) {
      DDs_rates_[3 * I + slip_direction] = slip_rate[I];
    }
  }
  void setSlipRateAlongSlipDirection(il::int_t &i, const double &slip_rate,
                                     il::int_t &slip_direction) {
    DDs_rates_[3 * i + slip_direction] = slip_rate;
  }
  void setTotalTractions(il::Array<double> &tot_tractions) {
    total_tractions_ = tot_tractions;
  }

  void setTotalTractions(il::int_t &elm_i, double &shear_stress_1,
                         double &shear_stress_2, double &normal_stress) {
    total_tractions_[3 * elm_i] = shear_stress_1;
    total_tractions_[3 * elm_i + 1] = shear_stress_2;
    total_tractions_[3 * elm_i + 2] = normal_stress;
  }
  void setCurrentTime(const double &curr_time) { current_time_ = curr_time; }
  void setTimeStep(const double &time_step) { time_step_ = time_step; }
  void setPressure(const il::Array<double> &pressure) { pressure_ = pressure; }
  void setPressure(il::int_t &i, const double &value) { pressure_[i] = value; }
  void setInitialPressure(const il::Array<double> &init_pressure) {
    pressure_o_ = init_pressure;
  }
};

/// RK45 Solution -> subclass of Solution class
class SolutionRK45 : public Solution {
 protected:
  il::Array<double> state_variables_;
  il::Array<double> plastic_fault_porosity_;
  il::Array<double> fault_hydraulic_aperture_;
  EQSim::Injection injection_;
  EQSim::FluidProperties fluid_properties_;
  EQSim::FaultProperties fault_properties_;
  EQSim::SolidMatrixProperties solidMatrix_Properties_;
  FrictionProperties *fric_properties_;
  PermeabilityProperties *permeab_properties_;
  std::string results_path_, baseFileName_;
  il::Array2D<double> ElastMatrix_;
  il::Array2D<il::int_t> neigh_elts_;

 public:
  SolutionRK45() = default;
  SolutionRK45(double &current_time, double &time_step,
               EQSim::SolverParameters &solver_parameters,
               EQSim::FluidProperties &fluid_prop,
               EQSim::SolidMatrixProperties &solidMatrix_Properties,
               il::Array<double> &insitu_tractions,
               il::Array<double> &total_tractions, il::Array<double> &DDs,
               il::Array<double> &DDs_rates, il::Array<double> &state_variables,
               il::Array<double> &plastic_fault_porosity,
               il::Array<double> &fault_hydraulic_aperture,
               il::Array<double> &pressure_o, il::Array<double> &pressure,
               EQSim::Injection &InjectionObj,
               EQSim::FaultProperties &fault_properties, EQSim::Mesh &mesh,
               FrictionProperties *fric_prop,
               PermeabilityProperties *permeab_prop, std::string &res_path,
               std::string &baseFileName, il::Array2D<double> const &elast_matrix,
               il::Array2D<il::int_t> &neigh_elts) {
    mesh_ = mesh;
    current_time_ = current_time;
    time_step_ = time_step;
    solver_parameters_ = solver_parameters;
    fluid_properties_ = fluid_prop;
    solidMatrix_Properties_ = solidMatrix_Properties;
    insitu_tractions_ = insitu_tractions;
    total_tractions_ = total_tractions;
    DDs_ = DDs;
    DDs_rates_ = DDs_rates;
    state_variables_ = state_variables;
    plastic_fault_porosity_ = plastic_fault_porosity;
    fault_hydraulic_aperture_ = fault_hydraulic_aperture;
    pressure_o_ = pressure_o;
    pressure_ = pressure;
    injection_ = InjectionObj;
    fault_properties_ = fault_properties;
    fric_properties_ = fric_prop;
    permeab_properties_ = permeab_prop;
    results_path_ = res_path;
    baseFileName_ = baseFileName;
    ElastMatrix_ = elast_matrix;
    neigh_elts_ = neigh_elts;
  }

  // Methods
  json createJsonObjectToExport(bool &export_background_stresses);
  void exportToJson(json &j_obj_to_export, std::string &filename);
  void exportToUBJson(json &j_obj_to_export, std::string &filename);

  // Getter methods
  Injection getInjectionObj() { return injection_; };
  FluidProperties getFluidProperties() { return fluid_properties_; };
  SolidMatrixProperties getSolidMatrixProperties() {
    return solidMatrix_Properties_;
  };
  FaultProperties getFaultProperties() { return fault_properties_; };
  FrictionProperties *getFrictionPtr() { return fric_properties_; };
  PermeabilityProperties *getPermeabPtr() { return permeab_properties_; };
  std::string getResPath() { return results_path_; };
  std::string getBaseFileName() { return baseFileName_; };
  il::Array2D<double> getElastMatrix() { return ElastMatrix_; };
  il::Array<double> getStateVariables() { return state_variables_; };
  il::Array<double> getPlasticFaultPorosity() {
    return plastic_fault_porosity_;
  };
  il::Array<double> getFaultHydraulicAperture() {
    return fault_hydraulic_aperture_;
  };
  double getInsituTractions(il::int_t i) const { return insitu_tractions_[i]; };
  il::Array2D<il::int_t> getNeighElts() const { return neigh_elts_; };

  // Setter methods
  void setStateVariables(il::Array<double> &state_variables) {
    state_variables_ = state_variables;
  };
  void setPlasticFaultPorosity(il::Array<double> &plastic_fault_porosity) {
    plastic_fault_porosity_ = plastic_fault_porosity;
  };
  void setFaultHydraulicAperture(il::Array<double> &fault_hydraulic_aperture) {
    fault_hydraulic_aperture_ = fault_hydraulic_aperture;
  };
};

}  // namespace EQSim
#endif  // INC_3DEQSIM_SRC_SOLUTION_H
