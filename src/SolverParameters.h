//
// Created by Federico Ciardo on 03.10.21. All rights reserved.
//

#ifndef INC_3DEQSIM_SRC_SOLVERPARAMETERS_H
#define INC_3DEQSIM_SRC_SOLVERPARAMETERS_H

#include <string>

//
// This structure collects global parameters controlling the numerical solver.
//

namespace EQSim {

struct SolverParameters {
  double initial_time;
  double time_Step;
  double maximum_time;
  double tolerance_RK = 10e-6;  // default value
  double time_step_amplification_factor;
  double time_step_reduction_factor;
  int export_current_solution_every_i_time_steps_;
  int export_stress_every_i_time_steps_;
  double loading_rate; //default value
  double adaptive_dt_constant{0.5}; // default value, used to adapt the time step based on the slip rate
  // Boundary condition applied to the fluid pressure field. Supported values
  // are "Dirichlet" and "Neumann". Default is Dirichlet to maintain backward
  // compatibility with previous input files.
  std::string pressure_boundary_condition{"Dirichlet"};
};

}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_SOLVERPARAMETERS_H
