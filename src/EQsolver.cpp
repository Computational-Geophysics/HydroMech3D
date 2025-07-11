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
                          double &Q, double &loading_rate);

///////////////////////////////////////////////////////////////////////////////

//                     EARTHQUAKE SOLVER

///////////////////////////////////////////////////////////////////////////////
void EQsolver(json &js, bool &checkRestart) {
  ///////////// IMPORT INPUT DATA FROM JSON CONFIG FILE /////////////

  std::string solver_description, basefilename, date, program_name, res_path;
  json j_injection, j_mesh, j_obsmesh, j_model_params, j_solver_params,
      j_fluid_params, j_fault_insitu_params,
      j_initial_conditions, j_rock_properties, j_fault_properties, j_fluid_flow,
      j_friction_properties, j_permeability_properties;

  /// IMPORT: #1 layer keystrings
  EQSim::ImportFirstLayerJsonInput(js, checkRestart, solver_description,
                                   basefilename, date, program_name, res_path,
                                   j_mesh, j_model_params, j_solver_params);

  /// IMPORT: #2 layer keystrings
  EQSim::ImportSecondLayerJsonInput(
      j_model_params, checkRestart, j_injection, j_fluid_params,
      j_fault_properties, j_fault_insitu_params,
      j_initial_conditions, j_rock_properties, j_fluid_flow);

  /// IMPORT: #3 layer keystrings
  EQSim::ImportThirdLayerJsonInput(
      checkRestart, j_fault_properties, j_friction_properties,
      j_permeability_properties);
  std::cout<<"All three layers are imported"<<std::endl;

  EQSim::FluidProperties FluidProperties = LoadFluidProperties(j_fluid_params);

  EQSim::SolidMatrixProperties SolidMatrixProperties = LoadMatrixProperties(j_rock_properties);

  EQSim::Mesh Mesh = LoadMeshData(j_mesh);

  EQSim::FaultProperties FaultProperties = LoadFaultProperties(
      Mesh, j_fault_properties, j_friction_properties,
      j_permeability_properties, j_initial_conditions);
  
  EQSim::FaultInSituStress FaultInSituStresses =
      LoadFaultInSituStressComponents(j_fault_insitu_params, Mesh);

  arma::vec insitu_tractions =
      FaultInSituStresses.AllInSituTractions(Mesh);
  std::cout<<"Insitu stress loaded"<<std::endl;

  EQSim::Injection Injection(j_injection);
  std::cout<<"Injection loaded"<<std::endl;

  EQSim::FrictionProperties *fric_coeff_properties =
      FaultProperties.getFrictionProperties();
  std::cout<<"Friction properties loaded"<<std::endl;

  EQSim::PermeabilityProperties *permeability_properties =
      FaultProperties.getPermeabilityProperties();
  std::cout<<"Permeability properties loaded"<<std::endl;

  EQSim::SolverParameters solver_parameters =
      LoadSolverParameters(j_solver_params);
  std::cout<<"Solver parameters loaded"<<std::endl;
  
  arma::wall_clock timer;
  timer.tic();
  // arma::mat ElastMatrix = EQSim::AssembleElastMat(Mesh, SolidMatrixProperties);
  arma::mat ElastMatrix{Mesh.getNumberOfDofs(), Mesh.getNumberOfDofs(), arma::fill::zeros};
  std::cout<<"Elastic matrix assembled"<<std::endl;
  double TAssemb = timer.toc();
  std::cout << "Time used to assemble elatic matrices in seconds: " << TAssemb << std::endl;

  arma::mat centroids = Mesh.getCentroids();
  arma::imat neigh_elts =
      Mesh.getNeighbourElements_UniformMesh(Mesh);
      // Mesh.getNeighbourElements_NonUniformMesh();
  std::cout<<"Initializing started"<<std::endl;

  // Initialization
  double current_time = solver_parameters.initial_time;
  double time_Step = solver_parameters.time_Step;
  double loading_rate = solver_parameters.loading_rate;
  arma::vec initial_DDs = FaultProperties.getInitialDDs();
  arma::vec initial_DDs_rates = FaultProperties.getInitialDDsRates();
  arma::vec initial_state_variables =
      FaultProperties.getInitialStateVariables();
  arma::vec initial_plastic_fault_porosity(Mesh.getNumberOfElts(), arma::fill::zeros);
  for (arma::uword I = 0; I < initial_plastic_fault_porosity.n_elem; ++I) {
    initial_plastic_fault_porosity[I] =
        FaultProperties.getInitialPlasticFaultPorosity(0);
  }
  arma::vec initial_fault_hydraulic_aperture(Mesh.getNumberOfElts(),
                                                     arma::fill::zeros);

  for (arma::uword I = 0; I < initial_fault_hydraulic_aperture.n_elem; ++I) {
    initial_fault_hydraulic_aperture[I] =
        FaultProperties.getInitialFaultHydraulicAperture(0);
  }
  arma::vec ambient_pressure =
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
  std::cout<<"RightHandSideODEs prepared"<<std::endl;

  EQSim::RK45 rk45(&RightHandSides, SolutionObj);
  std::cout<<"Initialization complete start to solve"<<std::endl;
  rk45.Solve();
  std::cout<<"RK45 solved"<<std::endl;
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
  double omega = (std::abs(slip_rate) * theta) / dc;
  //return -1. * omega * log(omega); // Slip law
  return 1 - omega; // Aging law
}

long double slip_rate_rhs(double &ref_slip_velocity, double &a_value,
                          double &b_value, double &theta, double &dc,
                          double &ref_fric_coeff, double &sigma_n, double &p,
                          double &slip_rate, double &GtimesSlipRate,
                          double &permeab, double &fluid_compress,
                          double &viscosity, double &hydraulic_aperture,
                          double &void_compress, double &porosity, double &d2_p,
                          double &Q, double &loading_rate) {
  // slip_rate = std::max(slip_rate, 1e-12);
  double shear_modulus = 3.0e10;    
  double c_s = 3000;                
  double eta = shear_modulus / ((1-0.25) * c_s); 
  double beta = porosity * (fluid_compress + void_compress);
  double omega = (std::abs(slip_rate) * theta) / dc;
  // double thetadot = -1. * omega * log(omega);
  double thetadot = 1 - omega;
  double pdot = ((permeab / (viscosity * beta)) * d2_p) +
                (Q / (hydraulic_aperture * beta));

  long double fric_coeff =
      ref_fric_coeff + (a_value * log(std::abs(slip_rate) / ref_slip_velocity)) +
      (b_value * log((theta * ref_slip_velocity) / dc));

  long double Numerator =
      std::abs(slip_rate) * ((theta * (fric_coeff * pdot - GtimesSlipRate - eta*slip_rate)) +
                            (b_value * (p - sigma_n) * thetadot));

  long double Denominator = ((a_value * sigma_n) - (a_value * p)) * theta;

  // std::cout<<"The value of Numerator is: "<<Numerator<<"The value of the Denominator is: "<<Denominator<<std::endl;
  // std::cout << "The value of GtimesSlipRate: "<<GtimesSlipRate<<"The value of the radiation damping is: "<<eta*slip_rate<<std::endl;
  if (std::isnan(Numerator)) {
    std::cerr << "Error: Numerator is NaN." << std::endl;
    std::cerr << "slip_rate = " << slip_rate << std::endl;
    std::cerr << "theta = " << theta << std::endl;
    std::cerr << "fric_coeff = " << fric_coeff << std::endl;
    std::cerr << "pdot = " << pdot << std::endl;
    std::cerr << "GtimesSlipRate = " << GtimesSlipRate << std::endl;
    std::cerr << "p = " << p << ", sigma_n = " << sigma_n << std::endl;
    std::cerr << "thetadot = " << thetadot << std::endl;
    assert(!std::isnan(Numerator));
}
  if (std::isnan(Denominator)) {
    std::cerr << "Error: Denominator is NaN." << std::endl;
    std::cerr << "a_value = " << a_value << ", sigma_n = " << sigma_n
              << ", p = " << p << ", theta = " << theta << std::endl;
    assert(!std::isnan(Denominator));
}
  long double dSlipRate_dt = Numerator / Denominator;
  // dSlipRate_dt += loading_rate;
  return dSlipRate_dt;
}

}  // namespace EQSim
