//
// Created by Federico Ciardo on 28.07.21. All rights reserved.
//

// Import from the project
#include "EQsolver.h"

#include "FrictionProperties.h"
#include "ImportJsonInputData.h"
#include "Injection.h"
#include "PermeabilityProperties.h"
#include "RK45.h"
#include "RightHandSidesODEs.h"
#include "Solution.h"
#include "AssembleElasticityMatrix.h"

namespace EQSim {

double pressure_rhs(double &permeab, double &fluid_compress,
                    double &viscosity, double &hydraulic_aperture,
                    double &void_compress, double &porosity, double &d2_p,
                    double &Q);

double state_rhs(double &theta, double &dc, double &slip_rate);

long double slip_rate_rhs(double &ref_slip_velocity, double &a_value,
                          double &b_value, double &theta, double &dc,
                          double &ref_fric_coeff, double &sigma_n, double &p,
                          double &slip_rate, double &GtimesSlipRate,
                          double &permeab, double &fluid_compress,
                          double &viscosity, double &hydraulic_aperture,
                          double &void_compress, double &porosity, double &d2_p,
                          double &Q);

///////////////////////////////////////////////////////////////////////////////

//                     EARTHQUAKE SOLVER

///////////////////////////////////////////////////////////////////////////////
void EQsolver(json &js, bool &checkRestart, il::io_t) {
  ///////////// IMPORT INPUT DATA FROM JSON CONFIG FILE /////////////

  std::string solver_description, basefilename, date, program_name, res_path;
  json j_injection, j_mesh, j_obsmesh, j_model_params, j_solver_params,
      j_fluid_params, j_fault_insitu_params,
      j_initial_conditions, j_rock_properties, j_fault_properties, j_fluid_flow,
      j_friction_properties, j_permeability_properties;

  /// IMPORT: #1 layer keystrings
  EQSim::ImportFirstLayerJsonInput(js, checkRestart, il::io, solver_description,
                                   basefilename, date, program_name, res_path,
                                   j_mesh, j_model_params, j_solver_params);

  /// IMPORT: #2 layer keystrings
  EQSim::ImportSecondLayerJsonInput(
      j_model_params, checkRestart, il::io, j_injection, j_fluid_params,
      j_fault_properties, j_fault_insitu_params,
      j_initial_conditions, j_rock_properties, j_fluid_flow);

  /// IMPORT: #3 layer keystrings
  EQSim::ImportThirdLayerJsonInput(
      checkRestart, j_fault_properties, il::io, j_friction_properties,
      j_permeability_properties);

  EQSim::FluidProperties FluidProperties = LoadFluidProperties(j_fluid_params);

  EQSim::SolidMatrixProperties SolidMatrixProperties = LoadMatrixProperties(j_rock_properties);

  EQSim::Mesh Mesh = LoadMeshData(j_mesh);

  EQSim::FaultProperties FaultProperties = LoadFaultProperties(
      Mesh, j_fault_properties, j_friction_properties,
      j_permeability_properties, j_initial_conditions);

  EQSim::FaultInSituStress FaultInSituStresses =
      LoadFaultInSituStressComponents(j_fault_insitu_params, Mesh);

  il::Array<double> insitu_tractions =
      FaultInSituStresses.AllInSituTractions(Mesh);

  EQSim::Injection Injection(j_injection);

  EQSim::FrictionProperties *fric_coeff_properties =
      FaultProperties.getFrictionProperties();
  EQSim::PermeabilityProperties *permeability_properties =
      FaultProperties.getPermeabilityProperties();

  EQSim::SolverParameters solver_parameters =
      LoadSolverParameters(j_solver_params);

  il::Array2D<double> ElastMatrix = EQSim::AssembleElastMat(Mesh, SolidMatrixProperties);

  il::Array2D<double> centroids = Mesh.getCentroids();
  il::Array2D<il::int_t> neigh_elts =
      Mesh.getNeighbourElements_UniformMesh(Mesh);

  // Initialization
  double current_time = solver_parameters.initial_time;
  double time_Step = solver_parameters.time_Step;

  il::Array<double> initial_DDs = FaultProperties.getInitialDDs();
  il::Array<double> initial_DDs_rates = FaultProperties.getInitialDDsRates();
  il::Array<double> initial_state_variables =
      FaultProperties.getInitialStateVariables();
  il::Array<double> initial_plastic_fault_porosity{Mesh.getNumberOfElts(), 0.};
  for (il::int_t I = 0; I < initial_plastic_fault_porosity.size(); ++I) {
    initial_plastic_fault_porosity[I] =
        FaultProperties.getInitialPlasticFaultPorosity(0);
  }
  il::Array<double> initial_fault_hydraulic_aperture{Mesh.getNumberOfElts(),
                                                     0.};

  for (il::int_t I = 0; I < initial_fault_hydraulic_aperture.size(); ++I) {
    initial_fault_hydraulic_aperture[I] =
        FaultProperties.getInitialFaultHydraulicAperture(0);
  }
  il::Array<double> ambient_pressure =
      FaultProperties.getAmbientPressureDistribution();

  EQSim::SolutionRK45 SolutionObj(
      current_time, time_Step, solver_parameters, FluidProperties,
      SolidMatrixProperties, insitu_tractions, insitu_tractions, initial_DDs,
      initial_DDs_rates, initial_state_variables,
      initial_plastic_fault_porosity, initial_fault_hydraulic_aperture,
      ambient_pressure, ambient_pressure, Injection, FaultProperties, Mesh,
      fric_coeff_properties, permeability_properties,
      res_path, basefilename, ElastMatrix, neigh_elts);

  EQSim::RightHandSideODEs RightHandSides(pressure_rhs, state_rhs,
                                          slip_rate_rhs);

  EQSim::RK45 rk45(&RightHandSides, SolutionObj);
  rk45.Solve();
}

double pressure_rhs(double &permeab, double &fluid_compress,
                    double &viscosity, double &hydraulic_aperture,
                    double &void_compress, double &porosity, double &d2_p,
                    double &Q) {
  double beta = porosity * (fluid_compress + void_compress);

  return ((permeab / (viscosity * beta)) * d2_p) +
         (Q / (hydraulic_aperture * beta));
}

double state_rhs(double &theta, double &dc, double &slip_rate) {
  double omega = (il::abs(slip_rate) * theta) / dc;

  return -1. * omega * log(omega);
}

long double slip_rate_rhs(double &ref_slip_velocity, double &a_value,
                          double &b_value, double &theta, double &dc,
                          double &ref_fric_coeff, double &sigma_n, double &p,
                          double &slip_rate, double &GtimesSlipRate,
                          double &permeab, double &fluid_compress,
                          double &viscosity, double &hydraulic_aperture,
                          double &void_compress, double &porosity, double &d2_p,
                          double &Q) {

  double beta = porosity * (fluid_compress + void_compress);
  double omega = (il::abs(slip_rate) * theta) / dc;
  double thetadot = -1. * omega * log(omega);
  double pdot = ((permeab / (viscosity * beta)) * d2_p) +
                (Q / (hydraulic_aperture * beta));

  long double fric_coeff =
      ref_fric_coeff + (a_value * log(il::abs(slip_rate) / ref_slip_velocity)) +
      (b_value * log((theta * ref_slip_velocity) / dc));

  long double Numerator =
      il::abs(slip_rate) * ((theta * (fric_coeff * pdot - GtimesSlipRate)) +
                            (b_value * (p - sigma_n) * thetadot));

  long double Denominator = ((a_value * sigma_n) - (a_value * p)) * theta;

  return Numerator / Denominator;
}

}  // namespace EQSim