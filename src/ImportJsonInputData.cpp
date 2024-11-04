//
// Created by Federico Ciardo on 28.07.21. All rights reserved.
//

#include "ImportJsonInputData.h"

#include "FrictionProperties.h"
#include "PermeabilityProperties.h"

namespace EQSim {

void ImportFirstLayerJsonInput(json &js, bool &checkRestart, il::io_t,
                               std::string &solver_description,
                               std::string &basefilename, std::string &date,
                               std::string &program_name, std::string &res_path,
                               json &j_mesh, json &j_model_params,
                               json &j_solver_params) {
  if (!checkRestart) {
    /// New analysis
    if (js.count("Date") == 1) {
      date = js["Date"].get<std::string>();
      std::cout << date << std::endl;
    } else {
      std::cout << "No date in input file - please add it" << std::endl;
      il::abort();
    }

    if (js.count("Program name") == 1) {
      program_name = js["Program name"].get<std::string>();
    } else {
      std::cout << "No program name in input file - please add it" << std::endl;
      il::abort();
    }

    if (js.count("Results path") == 1) {
      res_path = js["Results path"].get<std::string>();
    } else {
      std::cout << "No results path in input file - please add it" << std::endl;
      il::abort();
    }

    if (js.count("Base name") == 1) {
      basefilename = js["Base name"].get<std::string>();
    } else {
      std::cout << "No base name for output file(s)" << std::endl;
      il::abort();
    }

    if (js.count("Description") == 1) {
      solver_description = js["Description"].get<std::string>();
      std::cout << "-------------------------------------" << std::endl;
      std::cout << "3DEQSim: a " << solver_description << std::endl;
      std::cout << "-------------------------------------" << std::endl;
    } else {
      std::cout << "No solver description in input file - please add it"
                << std::endl;
      il::abort();
    }

    if (js.count("Mesh") == 1) {
      j_mesh = js["Mesh"];
    } else {
      std::cout << "No Mesh in json input file " << std::endl;
      il::abort();
    }

    if (js.count("Model parameters") == 1) {
      j_model_params = js["Model parameters"];
    } else {
      std::cout << "No model parameters in json input file " << std::endl;
      il::abort();
    }

    if (js.count("Solver parameters") == 1) {
      j_solver_params = js["Solver parameters"];
    } else {
      std::cout << "No solver parameters in json input file " << std::endl;
      il::abort();
    }

  } else {
    /// Restart analysis
    // TODO
  }
}

void ImportSecondLayerJsonInput(json &j_model_params, bool &checkRestart,
                                il::io_t, json &j_injection,
                                json &j_fluid_params, json &j_fault_properties,
                                json &j_fault_insitu_params,
                                json &j_initial_conditions,
                                json &j_rock_properties, json &j_fluid_flow) {
  if (!checkRestart) {
    /// New analysis

    if ((j_model_params.count("Fluid properties") == 1)) {
      j_fluid_params = j_model_params["Fluid properties"];
    } else {
      std::cout << "No fluid properties in json input file ";
      il::abort();
    }

    if ((j_model_params.count("Fault in-situ conditions") == 1)) {
      j_fault_insitu_params = j_model_params["Fault in-situ conditions"];
    } else {
      std::cout << "No fault in-situ conditions in json input file ";
      il::abort();
    }

    if ((j_model_params.count("Initial conditions") == 1)) {
      j_initial_conditions = j_model_params["Initial conditions"];
    } else {
      std::cout << "No initial conditions in json input file ";
      il::abort();
    }

    if ((j_model_params.count("Rock properties") == 1)) {
      j_rock_properties = j_model_params["Rock properties"];
    } else {
      std::cout << "No rock properties in json input file ";
      il::abort();
    }

    if ((j_model_params.count("Fault properties") == 1)) {
      j_fault_properties = j_model_params["Fault properties"];
    } else {
      std::cout << "No fault properties in json input file ";
      il::abort();
    }

    if ((j_model_params.count("Injection condition") == 1)) {
      j_injection = j_model_params["Injection condition"];
    } else {
      std::cout << "No injection condition in json input file ";
      il::abort();
    }

    if ((j_model_params.count("Fluid flow") == 1)) {
      j_fluid_flow = j_model_params["Fluid flow"];
    } else {
      std::cout << "No fluid flow in json input file ";
      il::abort();
    }

  } else {
    /// Restart analysis
    // TODO
  }
}

void ImportThirdLayerJsonInput(bool &checkRestart, json &j_fault_properties,
                               il::io_t, json &j_friction_properties,
                               json &j_permeability_properties) {
  if (!checkRestart) {
    /// New analysis

    if ((j_fault_properties.count("Friction") == 1)) {
      j_friction_properties = j_fault_properties["Friction"];
    } else {
      std::cout << "No friction properties in json input file ";
      il::abort();
    }

    if ((j_fault_properties.count("Permeability") == 1)) {
      j_permeability_properties = j_fault_properties["Permeability"];
    } else {
      std::cout << "No permeability properties in json input file ";
      il::abort();
    }

  } else {
    /// Restart analysis
    // TODO
  }
}

// Load FluidProperties Properties
EQSim::FluidProperties LoadFluidProperties(json &j_fluid_params) {
  std::string fluid_rheology;
  double fluid_density, fluid_compressibility, fluid_viscosity;

  if (j_fluid_params.count("Rheology") == 1) {
    fluid_rheology = j_fluid_params["Rheology"].get<std::string>();
  } else {
    std::cout << "No fluid rheology in input data!";
    il::abort();
  }

  if (j_fluid_params.count("Density") == 1) {
    fluid_density = j_fluid_params["Density"].get<double>();
  } else {
    std::cout << "No fluid density in input data!";
    il::abort();
  }

  if (j_fluid_params.count("Fluid Compressibility") == 1) {
    fluid_compressibility =
        j_fluid_params["Fluid Compressibility"].get<double>();
  } else {
    std::cout << "No Fluid Compressibility in input data!";
    il::abort();
  }

  if (j_fluid_params.count("Viscosity") == 1) {
    fluid_viscosity = j_fluid_params["Viscosity"].get<double>();
  } else {
    std::cout << "No fluid viscosity in input data!";
    il::abort();
  }

  // Call the constructor
  EQSim::FluidProperties Fluid(fluid_density, fluid_compressibility,
                               fluid_viscosity);
  return Fluid;
}

// Load Fault Properties
EQSim::FaultProperties LoadFaultProperties(EQSim::Mesh &Mesh,
                                           json &j_fault_properties,
                                           json &j_friction_properties,
                                           json &j_permeability_properties,
                                           json &j_initial_conditions) {
  double porosity, void_compressibility;

  if (j_fault_properties.count("Porosity") == 1) {
    porosity = j_fault_properties["Porosity"].get<double>();
  } else {
    std::cout << "No fault porosity in input data!";
    il::abort();
  }

  il::Array<double> initial_plastic_fault_porosity;
  if (j_fault_properties.count("Initial plastic porosity") == 1) {
    il::int_t NumberOFPlasticFaultPorosity =
        j_fault_properties["Initial plastic porosity"].size();
    initial_plastic_fault_porosity.Resize(NumberOFPlasticFaultPorosity);
    for (il::int_t I = 0; I < initial_plastic_fault_porosity.size(); ++I) {
      initial_plastic_fault_porosity[I] =
          j_fault_properties["Initial plastic porosity"][I];
    }
  } else {
    std::cout << "No Initial plastic fault porosity in input data!";
    il::abort();
  }

  il::Array<double> initial_hydraulic_aperture;
  if (j_fault_properties.count("Initial hydraulic aperture") == 1) {
    il::int_t NumberOFFaultHydraulicApertures =
        j_fault_properties["Initial hydraulic aperture"].size();
    initial_hydraulic_aperture.Resize(NumberOFFaultHydraulicApertures);
    for (il::int_t I = 0; I < initial_hydraulic_aperture.size(); ++I) {
      initial_hydraulic_aperture[I] =
          j_fault_properties["Initial hydraulic aperture"][I];
    }
  } else {
    std::cout << "No Initial hydraulic aperture in input data!";
    il::abort();
  }

  if (j_fault_properties.count("Void compressibility") == 1) {
    void_compressibility =
        j_fault_properties["Void compressibility"].get<double>();
  } else {
    std::cout << "No void compressibility in input data!";
    il::abort();
  }

  EQSim::FrictionProperties *fric_coeff_properties =
      nullptr;  // Initialization with null pointer

  if (j_friction_properties.count("Type") == 1) {
    if (j_friction_properties["Type"].get<std::string>() ==
        "Rate and State Friction") {
      fric_coeff_properties = new RateAndStateFriction;
    } else {
      std::cout << "Type of friction coefficient is not valid! " << std::endl;
      il::abort();
    }
  } else {
    std::cout << "Type of friction coefficient is not specified in input file! "
              << std::endl;
    il::abort();
  }

  // Now we declare a pointer for permeability that point to the type of
  // permeability law the user want to use. Such a specification must be
  // inserted in the input file, else the simulation is aborted.
  EQSim::PermeabilityProperties *permeability_properties =
      nullptr;  // Initialization with null pointer

  if (j_permeability_properties.count("Type") == 1) {
    if (j_permeability_properties["Type"].get<std::string>() == "Constant") {
      permeability_properties = new ConstantPermeability;
    }
    else {
      std::cout << "Type of permeability law is not valid! " << std::endl;
      il::abort();
    }
  } else {
    std::cout << "Type of permeability is not specified in input file! "
              << std::endl;
    il::abort();
  }

  // Now that the pointer points to the desired law, set the parameters
  fric_coeff_properties->setFrictionParameters(j_friction_properties, Mesh);
  permeability_properties->setPermeabilityParameters(j_permeability_properties,
                                                     Mesh);

  // Import fault initial conditions in terms of DDs, DDs rate and state
  // variable
  if (j_initial_conditions.count("Initial DDs") != 1) {
    std::cout << "Initial DDs keyword is wrong in input file! " << std::endl;
    il::abort();
  }

  if (j_initial_conditions.count("Initial DDs rates") != 1) {
    std::cout << "Initial DDs rates keyword is wrong in input file! "
              << std::endl;
    il::abort();
  }

  if (j_initial_conditions.count("Initial state variable") != 1) {
    std::cout << "Initial state variable is wrong in input file! " << std::endl;
    il::abort();
  }

  if (j_initial_conditions.count("Internal loads") != 1) {
    std::cout << "Internal loads is wrong in input file! " << std::endl;
    il::abort();
  }

  IL_ASSERT(j_initial_conditions["Initial DDs"].size() ==
            Mesh.getNumberOfDofs());
  IL_ASSERT(j_initial_conditions["Initial DDs rates"].size() ==
            Mesh.getNumberOfDofs());
  IL_ASSERT(j_initial_conditions["Initial state variable"].size() ==
            Mesh.getNumberOfElts());
  IL_ASSERT(j_initial_conditions["Internal loads"].size() ==
            Mesh.getNumberOfDofs());
  il::Array<double> initial_DDs{Mesh.getNumberOfDofs(), 0.};
  il::Array<double> initial_DDs_rates{Mesh.getNumberOfDofs(), 0.};
  il::Array<double> initial_state_variables{Mesh.getNumberOfElts(), 0.};
  il::Array<double> internal_loads{Mesh.getNumberOfDofs(), 0.};
  for (il::int_t I = 0; I < Mesh.getNumberOfElts(); ++I) {
    initial_DDs[3 * I] = j_initial_conditions["Initial DDs"][3 * I];
    initial_DDs[3 * I + 1] = j_initial_conditions["Initial DDs"][3 * I + 1];
    initial_DDs[3 * I + 2] = j_initial_conditions["Initial DDs"][3 * I + 2];
    initial_DDs_rates[3 * I] = j_initial_conditions["Initial DDs rates"][3 * I];
    initial_DDs_rates[3 * I + 1] =
        j_initial_conditions["Initial DDs rates"][3 * I + 1];
    initial_DDs_rates[3 * I + 2] =
        j_initial_conditions["Initial DDs rates"][3 * I + 2];
    initial_state_variables[I] =
        j_initial_conditions["Initial state variable"][I];
    internal_loads[3 * I] = j_initial_conditions["Internal loads"][3 * I];
    internal_loads[3 * I + 1] =
        j_initial_conditions["Internal loads"][3 * I + 1];
    internal_loads[3 * I + 2] =
        j_initial_conditions["Internal loads"][3 * I + 2];
  }

  // Call the constructor
  EQSim::FaultProperties FaultProperties(
      porosity, void_compressibility, initial_hydraulic_aperture, fric_coeff_properties,
      permeability_properties, initial_DDs,
      initial_DDs_rates, initial_state_variables,
      initial_plastic_fault_porosity, internal_loads);

  return FaultProperties;
}

// Load Matrix Properties
EQSim::SolidMatrixProperties LoadMatrixProperties(json &j_rock_properties) {
  double YoungModulus, PoissonRatio, RockDensity;

  if (j_rock_properties.count("Young's modulus") == 1) {
    YoungModulus = j_rock_properties["Young's modulus"].get<double>();
  } else {
    std::cout << "No Young's modulus in input data!";
    il::abort();
  }

  if (j_rock_properties.count("Poisson's ratio") == 1) {
    PoissonRatio = j_rock_properties["Poisson's ratio"].get<double>();
  } else {
    std::cout << "No Poisson's ratio in input data!";
    il::abort();
  }

  if (j_rock_properties.count("Density") == 1) {
    RockDensity = j_rock_properties["Density"].get<double>();
  } else {
    std::cout << "No rock density in input data!";
    il::abort();
  }

  EQSim::SolidMatrixProperties MatrixProperties(YoungModulus, PoissonRatio,
                                                RockDensity);
  return MatrixProperties;
}

// Load Solver parameters
EQSim::SolverParameters LoadSolverParameters(json &j_solver_params) {
  EQSim::SolverParameters solver_parameter_struct;

  // Import some solver parameters
  if (j_solver_params.count("Initial time") == 1) {
    solver_parameter_struct.initial_time =
        j_solver_params["Initial time"].get<double>();
  } else {
    std::cout << "No Initial time in json input file ";
    il::abort();
  }

  if (j_solver_params.count("Time step") == 1) {
    solver_parameter_struct.time_Step =
        j_solver_params["Time step"].get<double>();
  } else {
    std::cout << "No Time step in json input file ";
    il::abort();
  }

  if (j_solver_params.count("Maximum time") == 1) {
    solver_parameter_struct.maximum_time =
        j_solver_params["Maximum time"].get<double>();
  } else {
    std::cout << "No Maximum in json input file ";
    il::abort();
  }

  if (j_solver_params.count("Time step amplification factor") == 1) {
    solver_parameter_struct.time_step_amplification_factor =
        j_solver_params["Time step amplification factor"].get<double>();
  } else {
    std::cout << "No Time step amplification factor in json input file ";
    il::abort();
  }

  if (j_solver_params.count("Time step reduction factor") == 1) {
    solver_parameter_struct.time_step_reduction_factor =
        j_solver_params["Time step reduction factor"].get<double>();
  } else {
    std::cout << "No Time step reduction factor in json input file ";
    il::abort();
  }

  if (j_solver_params.count("Tolerance RK") == 1) {
    solver_parameter_struct.tolerance_RK =
        j_solver_params["Tolerance RK"].get<double>();
  } else {
    std::cout << "No Tolerance RK in json input file - the default value will "
                 "be used";
  }

  if (j_solver_params.count("Export current solution every i time steps") ==
      1) {
    solver_parameter_struct.export_current_solution_every_i_time_steps_ =
        j_solver_params["Export current solution every i time steps"]
            .get<int>();
  } else {
    std::cout
        << "No Export current solution every i time steps in json input file ";
    il::abort();
  }

  return solver_parameter_struct;
}

// Load Mesh Data
EQSim::Mesh LoadMeshData(json &j_mesh) {
  if (j_mesh.count("Node coordinates") != 1) {
    std::cout << "No node coordinates in input data!";
    il::abort();
  }

  if (j_mesh.count("Connectivity Nodes") != 1) {
    std::cout << "No connectivity nodes in input data!";
    il::abort();
  }

  il::int_t Nelts = j_mesh["Connectivity Nodes"].size();
  il::int_t Nnodes = j_mesh["Node coordinates"].size();

  IL_EXPECT_FAST(Nelts > 0);
  IL_EXPECT_FAST(Nnodes > 0);
  IL_EXPECT_FAST(j_mesh["Node coordinates"][0].size() == 3);
  IL_EXPECT_FAST(j_mesh["Connectivity Nodes"][0].size() == 4);

  il::Array2D<double> node_coordinates{Nnodes, 3, 0.};
  il::Array2D<il::int_t> connettivity{Nelts, 4, 0};
  il::int_t interpolation_order = 0;

  for (il::int_t i = 0; i < Nnodes; i++) {
    node_coordinates(i, 0) = j_mesh["Node coordinates"][i][0];
    node_coordinates(i, 1) = j_mesh["Node coordinates"][i][1];
    node_coordinates(i, 2) = j_mesh["Node coordinates"][i][2];
  }

  for (il::int_t i = 0; i < Nelts; i++) {
    connettivity(i, 0) = j_mesh["Connectivity Nodes"][i][0];
    connettivity(i, 1) = j_mesh["Connectivity Nodes"][i][1];
    connettivity(i, 2) = j_mesh["Connectivity Nodes"][i][2];
    connettivity(i, 3) = j_mesh["Connectivity Nodes"][i][3];
  }

  if (j_mesh.count("Interpolation order") == 1) {
    interpolation_order = j_mesh["Interpolation order"].get<il::int_t>();
  } else {
    std::cout << "No interpolation order in input data!";
    il::abort();
  }

  // Call mesh constructor
  EQSim::Mesh Mesh(node_coordinates, connettivity,
                   interpolation_order);
  return Mesh;
}

EQSim::FaultInSituStress LoadFaultInSituStressComponents(json &j_fault_insitu,
                                                         EQSim::Mesh &Mesh) {
  std::string TypeInSitu;

  if (j_fault_insitu.count("Type") == 1) {
    TypeInSitu = j_fault_insitu["Type"].get<std::string>();
  } else {
    std::cout << "No type of in situ stress in input data!";
    il::abort();
  }

  il::int_t sxx_present = j_fault_insitu.count("Sxx");
  il::int_t syy_present = j_fault_insitu.count("Syy");
  il::int_t szz_present = j_fault_insitu.count("Szz");
  il::int_t sxy_present = j_fault_insitu.count("Sxy");
  il::int_t sxz_present = j_fault_insitu.count("Sxz");
  il::int_t syz_present = j_fault_insitu.count("Syz");

  if (sxx_present != 1 || syy_present != 1 || szz_present != 1 ||
      sxy_present != 1 || sxz_present != 1 || syz_present != 1) {
    std::cout << " error in input file - No fault stress values \n";
    il::abort();
  }

  il::int_t nxx = j_fault_insitu["Sxx"].size();
  il::int_t nyy = j_fault_insitu["Syy"].size();
  il::int_t nzz = j_fault_insitu["Szz"].size();
  il::int_t nxy = j_fault_insitu["Sxy"].size();
  il::int_t nxz = j_fault_insitu["Sxz"].size();
  il::int_t nyz = j_fault_insitu["Syz"].size();

  il::int_t Nelts = Mesh.getNumberOfElts();

  IL_EXPECT_FAST(Nelts == nxx && Nelts == nyy && Nelts == nzz && Nelts == nxy &&
                 Nelts == nxz && Nelts == nyz);

  il::Array<double> Sxx(Nelts, 0.), Syy(Nelts, 0.), Szz(Nelts, 0.),
      Sxy(Nelts, 0.), Sxz(Nelts, 0.), Syz(Nelts, 0.);

  for (il::int_t i = 0; i < Nelts; i++) {
    Sxx[i] = j_fault_insitu["Sxx"][i];
    Syy[i] = j_fault_insitu["Syy"][i];
    Szz[i] = j_fault_insitu["Szz"][i];
    Sxy[i] = j_fault_insitu["Sxy"][i];
    Sxz[i] = j_fault_insitu["Sxz"][i];
    Syz[i] = j_fault_insitu["Syz"][i];
  }

  FaultInSituStress FaultInSituStress(Sxx, Syy, Szz, Sxy, Sxz, Syz);
  return FaultInSituStress;
}

}  // namespace EQSim