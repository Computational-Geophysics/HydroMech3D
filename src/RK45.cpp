// Created by Federico Ciardo on 30.08.21. All rights reserved.

// Inclusion from Standard library
#include <iostream>
#include <armadillo> // Include Armadillo library
#include "FrictionProperties.h"
#include "PermeabilityProperties.h"
#include "RK45.h"
#include "SolidMatrixProperties.h"
#include "Utils.h"
#include "Hmatrix_bigwham.h"
#include "ExportBackgroundStress.h"
#include "ExportBackgroundStressDirect.h"
#include <limits>
namespace EQSim {
void RK45::Solve() {
    bool use_hmatrix = true;
    Mesh FaultMesh = Solution_.getMesh();
    arma::uword Nelts = FaultMesh.getNumberOfElts();
    arma::mat ElastMat = Solution_.getElastMatrix();
    // ElastMat.save("ElastMat.csv", arma::csv_ascii);
    // ElastMat.save("ElastMat.bin", arma::arma_binary);
    // ElastMat.save("ElastMat.h5", arma::hdf5_binary);
    double E  = Solution_.getSolidMatrixProperties().getYoungModulus();
    double nu = Solution_.getSolidMatrixProperties().getPoissonRatio();
    ExportBackgroundStressDirect stressExporter(FaultMesh, E, nu);

    // build once
    EQSim::HMatrixBigwham hmat(FaultMesh, E, nu);
    arma::mat coords(Nelts, 3);
    for (arma::uword i = 0; i < Nelts; ++i) {
        EQSim::ElementData element_data = FaultMesh.getElementData(i);  // Get node coordinates and centroid
        coords.row(i) = element_data.getCentroidElt().t();  
    }
    // Construct H-Matrix from ElastMat
    std::cout << "coords size = " 
              << coords.n_rows << " x " 
              << coords.n_cols << std::endl;   
    arma::uvec idx = arma::regspace<arma::uvec>(0, 3, ElastMat.n_rows - 1);


    EQSim::PermeabilityProperties *PermPtr = Solution_.getPermeabPtr();
    EQSim::FrictionProperties *FricPtr = Solution_.getFrictionPtr();
    EQSim::SolidMatrixProperties SolidMatrixProperties = Solution_.getSolidMatrixProperties();
    EQSim::FluidProperties FluidProperties = Solution_.getFluidProperties();
    double viscosity = FluidProperties.getFluidViscosity();
    double fluid_compressibility = FluidProperties.getFluidCompressibility();
    // Get solver parameters
    EQSim::SolverParameters solver_parameters = Solution_.getSolverParametersStructure();
    bool neumann_bc = (solver_parameters.pressure_boundary_condition == "Neumann");
    // Get Injection properties
    Injection InjectionObj = Solution_.getInjectionObj();
    // Get Fault properties & initial conditions from solution obj
    FaultProperties FaultProperties = Solution_.getFaultProperties();
    double fault_porosity = FaultProperties.getFaultPorosity();
    double void_compressibility = FaultProperties.getFaultVoidCompressibility();
    // Get neighbour 2D array & get the boundary elements
    arma::imat neigh_elts = Solution_.getNeighElts();
    // Get boundary elements
    arma::ivec boundary_elts = FaultMesh.getBoundaryElements();
    arma::uword slip_direction = 0;
    arma::uword Aux = 0;
    if (Solution_.getInsituTractions(Aux) == 0) {
        slip_direction = 1;
    }
    // Initialization
    /// Eq1 -> solve for pressure
    arma::vec average_slope1_eq1(Nelts, arma::fill::zeros);
    arma::vec average_slope2_eq1(Nelts, arma::fill::zeros);
    arma::vec average_slope3_eq1(Nelts, arma::fill::zeros);
    arma::vec average_slope4_eq1(Nelts, arma::fill::zeros);
    arma::vec average_slope5_eq1(Nelts, arma::fill::zeros);
    arma::vec average_slope6_eq1(Nelts, arma::fill::zeros);
    /// Eq2 -> solve for slip
    arma::vec average_slope1_eq2(Nelts, arma::fill::zeros);
    arma::vec average_slope2_eq2(Nelts, arma::fill::zeros);
    arma::vec average_slope3_eq2(Nelts, arma::fill::zeros);
    arma::vec average_slope4_eq2(Nelts, arma::fill::zeros);
    arma::vec average_slope5_eq2(Nelts, arma::fill::zeros);
    arma::vec average_slope6_eq2(Nelts, arma::fill::zeros);
    /// Eq3 -> solve for state variable
    arma::vec average_slope1_eq3(Nelts, arma::fill::zeros);
    arma::vec average_slope2_eq3(Nelts, arma::fill::zeros);
    arma::vec average_slope3_eq3(Nelts, arma::fill::zeros);
    arma::vec average_slope4_eq3(Nelts, arma::fill::zeros);
    arma::vec average_slope5_eq3(Nelts, arma::fill::zeros);
    arma::vec average_slope6_eq3(Nelts, arma::fill::zeros);
    /// Eq4 -> solve for slip rate
    arma::vec average_slope1_eq4(Nelts, arma::fill::zeros);
    arma::vec average_slope2_eq4(Nelts, arma::fill::zeros);
    arma::vec average_slope3_eq4(Nelts, arma::fill::zeros);
    arma::vec average_slope4_eq4(Nelts, arma::fill::zeros);
    arma::vec average_slope5_eq4(Nelts, arma::fill::zeros);
    arma::vec average_slope6_eq4(Nelts, arma::fill::zeros);
    arma::vec Q_vector;

    // Helper lambda to enforce pressure boundary conditions for a given
    // element index according to the solver settings.
    auto apply_pressure_bc = [&](arma::uword idx) {
        if (solver_parameters.pressure_boundary_condition == "Dirichlet") {
            Solution_.setPressure(idx, Solution_.getInitPressure(idx));
        } else if (solver_parameters.pressure_boundary_condition == "Neumann") {
            double sum = 0.0;
            arma::uword count = 0;
            for (arma::uword j = 0; j < neigh_elts.n_cols; ++j) {
                int n = neigh_elts(idx, j);
                if (n >= 0 && n != static_cast<int>(idx)) {
                    sum += Solution_.getPressure(n);
                    ++count;
                }
            }
            if (count > 0) {
                Solution_.setPressure(idx, sum / static_cast<double>(count));
            }
        }
    };

    auto apply_pressure_bc_vec = [&](arma::vec& vec, arma::uword idx) {
        if (solver_parameters.pressure_boundary_condition == "Dirichlet") {
            vec[idx] = Solution_.getInitPressure(idx);
        } else if (solver_parameters.pressure_boundary_condition == "Neumann") {
            double sum = 0.0;
            arma::uword count = 0;
            for (arma::uword j = 0; j < neigh_elts.n_cols; ++j) {
                int n = neigh_elts(idx, j);
                if (n >= 0 && n != static_cast<int>(idx)) {
                    sum += vec[n];
                    ++count;
                }
            }
            if (count > 0) {
                vec[idx] = sum / static_cast<double>(count);
            }
        }
    };
    double perm_I, dc_I, theta_I, slip_rate_I, sigma_n_I, ref_fric_I, a_I, b_I, ref_slip_rate_I, press_I, epsilon_dilat_I, ElastMatTimesDDsRates_I, plastic_porosity_I;
    double Q_i;
    double hydraulic_aperture_I;
    double loading_rate = solver_parameters.loading_rate;
    arma::vec initial_hydraulic_aperture = Solution_.getFaultHydraulicAperture();
    arma::vec ElastMatTimesDDs;
    arma::vec ElastMatTimesDDsRates;
    arma::vec laplacian_pressure;
    arma::vec previous_pressure;
    arma::vec previous_state_variables;
    arma::vec previous_plastic_fault_porosity;
    arma::vec previous_fault_hydraulic_aperture;
    arma::vec previous_DDs;
    arma::vec previous_DDs_Rates;
    arma::vec new_pressure_order4(Nelts, arma::fill::zeros);
    arma::vec new_pressure_order5(Nelts, arma::fill::zeros);
    arma::vec new_theta_order4(Nelts, arma::fill::zeros);
    arma::vec new_theta_order5(Nelts, arma::fill::zeros);
    arma::vec new_plastic_fault_porosity_order4(Nelts, arma::fill::zeros);
    arma::vec new_plastic_fault_porosity_order5(Nelts, arma::fill::zeros);
    arma::vec new_hydraulic_aperture_order4(Nelts, arma::fill::zeros);
    arma::vec new_hydraulic_aperture_order5(Nelts, arma::fill::zeros);
    arma::vec new_slip_order4(Nelts, arma::fill::zeros);
    arma::vec new_slip_order5(Nelts, arma::fill::zeros);
    arma::vec new_slip_rates_order4(Nelts, arma::fill::zeros);
    arma::vec new_slip_rates_order5(Nelts, arma::fill::zeros);
    arma::vec diff_pressure45(Nelts, arma::fill::zeros);
    arma::vec diff_slip_rates45(Nelts, arma::fill::zeros);
    arma::vec diff_theta45(Nelts, arma::fill::zeros);
    double error_pressure45, error_slip_rate45, error_theta45, max_error;

    arma::uword iter;
    std::string filenameSolution;
    std::string ncfilenameSolution;

    json json_export_solution;
    double time = Solution_.getCurrentTime();
    arma::uword counter = 0;
    arma::mat stresses_on_background_mesh;
    bool export_background_stresses = true;

    arma::vec Vload_vec(3*Nelts, arma::fill::zeros);
    Vload_vec.elem(idx).fill(solver_parameters.loading_rate);
    arma::vec hmatVload   = hmat.multiply(Vload_vec);
    double minTx = hmatVload.elem(idx).min();
    arma::vec constant_stress_loading(3*Nelts, arma::fill::zeros);
    constant_stress_loading.elem(idx).fill(2*minTx);
    std::cout << "Constant stress loading rate: " << minTx << "Pa/s" <<std::endl;
    
    arma::uword srcIdx = InjectionObj.getSourceElement1();   // l‑value
    std::cout << "***** Start time marching procedure *****" << std::endl;
    while (time <= solver_parameters.maximum_time) {
        ++counter;
        std::cout << "- Current time: " << time << " s" << std::endl;
        previous_pressure = Solution_.getPressure();
        previous_state_variables = Solution_.getStateVariables();
        previous_plastic_fault_porosity = Solution_.getPlasticFaultPorosity();
        previous_fault_hydraulic_aperture = Solution_.getFaultHydraulicAperture();
        previous_DDs = Solution_.getDDs();
        previous_DDs_Rates = Solution_.getDDsRates();
        Q_vector = InjectionObj.calculate_Q_vector(FaultMesh, FluidProperties, time);
        //Q_vector = InjectionObj.calculate_Q_vector(FaultMesh, FluidProperties, time, previous_pressure, Solution_.getTimeStep());
        error_pressure45 = 2.;
        error_slip_rate45 = 2.;
        error_theta45 = 2.;
        iter = 0;
        while ((error_pressure45 > solver_parameters.tolerance_RK) || (error_slip_rate45 > solver_parameters.tolerance_RK) || (error_theta45 > solver_parameters.tolerance_RK)) {
            ++iter;
            laplacian_pressure = EQSim::CalculatePressureLaplacianUniformMeshOnly(neigh_elts, FaultMesh, previous_pressure, neumann_bc);
            // ElastMatTimesDDsRates = ElastMat * previous_DDs_Rates;
            // Convert arma::vec to std::vector<double>
            // ElastMatTimesDDsRates = arma::conv_to<arma::vec>::from(
            //     hmatrix->multiply(previous_DDs_Rates)
            // );
            /* std::cout<<"Size for elastmat: "<<ElastMat.n_rows<<"x"<<ElastMat.n_cols<<std::endl;
            std::cout << "previousDDrates size = " 
                << previous_DDs_Rates.n_rows << " x " 
                << previous_DDs_Rates.n_cols << std::endl;*/   
            // ElastMatTimesDDsRates = hmatrix.multiply(previous_DDs_Rates); //Using Hmatrix
            if (use_hmatrix) {
                // Use the H-matrix multiplication
                // ElastMatTimesDDsRates = hmat.multiply(previous_DDs_Rates + Vload_vec);
                ElastMatTimesDDsRates = hmat.multiply(previous_DDs_Rates);
                ElastMatTimesDDsRates -= constant_stress_loading;
            } else {
                // Use the dense multiplication
                ElastMatTimesDDsRates = ElastMat * (previous_DDs_Rates - Vload_vec);
            }
            for (arma::uword I = 0; I < Nelts; ++I) {
                perm_I = PermPtr->getPermeabilities(I);
                dc_I = FricPtr->get_state_evolution_distances_(I);
                Q_i = Q_vector[I];
                hydraulic_aperture_I = previous_fault_hydraulic_aperture[I];
                press_I = previous_pressure[I];
                slip_rate_I = previous_DDs_Rates[3 * I + slip_direction];
                theta_I = previous_state_variables[I];
                plastic_porosity_I = previous_plastic_fault_porosity[I];
                ElastMatTimesDDsRates_I = ElastMatTimesDDsRates[3 * I + slip_direction];
                average_slope1_eq1[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc1_(perm_I, fluid_compressibility, viscosity, hydraulic_aperture_I,
                void_compressibility, fault_porosity, laplacian_pressure[I], Q_i));
                Solution_.setPressure(I, previous_pressure[I] + ((1. / 4.) * average_slope1_eq1[I]));
                //if (InjectionObj.isHoldingPhase()){
                //    Solution_.setPressure( srcIdx, InjectionObj.getTargetPressure());
                //    laplacian_pressure[srcIdx] = 0.0;}
                average_slope1_eq2[I] = Solution_.getTimeStep() * slip_rate_I;
                average_slope1_eq3[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc2_(theta_I, dc_I, slip_rate_I));
                sigma_n_I = Solution_.getNormalTotalTractions(I);
                ref_fric_I = FricPtr->get_reference_friction_coefficients_(I);
                a_I = FricPtr->get_a_values(I);
                b_I = FricPtr->get_b_values(I);
                ref_slip_rate_I = FricPtr->get_reference_velocities(I);
                average_slope1_eq4[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc3_(ref_slip_rate_I, a_I, b_I, theta_I, dc_I, ref_fric_I, sigma_n_I, press_I, slip_rate_I, ElastMatTimesDDsRates_I, perm_I, fluid_compressibility, viscosity, hydraulic_aperture_I, void_compressibility, fault_porosity, laplacian_pressure[I], Q_i, loading_rate));
                Solution_.setSlipRateAlongSlipDirection(I, previous_DDs_Rates[3 * I + slip_direction] + ((1. / 4.) * average_slope1_eq4[I]), slip_direction);
            }
            for (arma::uword be = 0; be < boundary_elts.n_elem; ++be) {
                apply_pressure_bc(boundary_elts[be]);
            }
            laplacian_pressure = EQSim::CalculatePressureLaplacianUniformMeshOnly(neigh_elts, FaultMesh, Solution_.getPressure(), neumann_bc);
            if (use_hmatrix) {
                // Use the H-matrix multiplication
                // ElastMatTimesDDsRates = hmat.multiply(Solution_.getDDsRates() + Vload_vec);
                ElastMatTimesDDsRates = hmat.multiply(Solution_.getDDsRates());
                ElastMatTimesDDsRates += constant_stress_loading;
            } else {
                ElastMatTimesDDsRates = ElastMat * (Solution_.getDDsRates() - Vload_vec);
            }
            for (arma::uword I = 0; I < Nelts; ++I) {
                perm_I = PermPtr->getPermeabilities(I);
                dc_I = FricPtr->get_state_evolution_distances_(I);
                Q_i = Q_vector[I];
                hydraulic_aperture_I = previous_fault_hydraulic_aperture[I];
                press_I = Solution_.getPressure(I);
                slip_rate_I = Solution_.getSlipRateAlongSlipDirection(I, slip_direction);
                theta_I = previous_state_variables[I] + ((1. / 4.) * average_slope1_eq3[I]);
                plastic_porosity_I = previous_plastic_fault_porosity[I];
                ElastMatTimesDDsRates_I = ElastMatTimesDDsRates[3 * I + slip_direction];
                average_slope2_eq1[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc1_(perm_I, fluid_compressibility, viscosity, hydraulic_aperture_I, void_compressibility, fault_porosity, laplacian_pressure[I], Q_i));
                Solution_.setPressure(I, previous_pressure[I] + ((3. / 32.) * average_slope1_eq1[I]) + ((9. / 32.) * average_slope2_eq1[I]));
                //if (InjectionObj.isHoldingPhase()){
                //    Solution_.setPressure( srcIdx, InjectionObj.getTargetPressure());
                //    laplacian_pressure[srcIdx] = 0.0;}
                average_slope2_eq2[I] = Solution_.getTimeStep() * slip_rate_I;
                average_slope2_eq3[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc2_(theta_I, dc_I, slip_rate_I));
                sigma_n_I = Solution_.getNormalTotalTractions(I);
                ref_fric_I = FricPtr->get_reference_friction_coefficients_(I);
                a_I = FricPtr->get_a_values(I);
                b_I = FricPtr->get_b_values(I);
                ref_slip_rate_I = FricPtr->get_reference_velocities(I);
                average_slope2_eq4[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc3_(ref_slip_rate_I, a_I, b_I, theta_I, dc_I, ref_fric_I, sigma_n_I, press_I, slip_rate_I, ElastMatTimesDDsRates_I, perm_I, fluid_compressibility, viscosity, hydraulic_aperture_I, void_compressibility, fault_porosity, laplacian_pressure[I], Q_i, loading_rate));
                Solution_.setSlipRateAlongSlipDirection(I, previous_DDs_Rates[3 * I + slip_direction] + ((3. / 32.) * average_slope1_eq4[I]) + ((9. / 32.) * average_slope2_eq4[I]), slip_direction);
            }
            for (arma::uword be = 0; be < boundary_elts.n_elem; ++be) {
                apply_pressure_bc(boundary_elts[be]);
            }
            laplacian_pressure = EQSim::CalculatePressureLaplacianUniformMeshOnly(neigh_elts, FaultMesh, Solution_.getPressure(), neumann_bc);
            if (use_hmatrix) {
                // Use the H-matrix multiplication
                // ElastMatTimesDDsRates = hmat.multiply(Solution_.getDDsRates() + Vload_vec);
                ElastMatTimesDDsRates = hmat.multiply(Solution_.getDDsRates());
                ElastMatTimesDDsRates -= constant_stress_loading;
            } else {
                ElastMatTimesDDsRates = ElastMat * (Solution_.getDDsRates() - Vload_vec);
            }
            for (arma::uword I = 0; I < Nelts; ++I) {
                perm_I = PermPtr->getPermeabilities(I);
                dc_I = FricPtr->get_state_evolution_distances_(I);
                Q_i = Q_vector[I];
                hydraulic_aperture_I = previous_fault_hydraulic_aperture[I];
                press_I = Solution_.getPressure(I);
                slip_rate_I = Solution_.getSlipRateAlongSlipDirection(I, slip_direction);
                theta_I = previous_state_variables[I] + ((3. / 32.) * average_slope1_eq3[I]) + ((9. / 32.) * average_slope2_eq3[I]);
                plastic_porosity_I = previous_plastic_fault_porosity[I];
                ElastMatTimesDDsRates_I = ElastMatTimesDDsRates[3 * I + slip_direction];
                average_slope3_eq1[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc1_(perm_I, fluid_compressibility, viscosity, hydraulic_aperture_I, void_compressibility, fault_porosity, laplacian_pressure[I], Q_i));
                Solution_.setPressure(I, previous_pressure[I] + ((1932. / 2197.) * average_slope1_eq1[I]) - ((7200. / 2197.) * average_slope2_eq1[I]) + ((7296. / 2197.) * average_slope3_eq1[I]));
                //if (InjectionObj.isHoldingPhase()){
                //    Solution_.setPressure( srcIdx, InjectionObj.getTargetPressure());
                //    laplacian_pressure[srcIdx] = 0.0;}
                average_slope3_eq2[I] = Solution_.getTimeStep() * slip_rate_I;
                average_slope3_eq3[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc2_(theta_I, dc_I, slip_rate_I));
                sigma_n_I = Solution_.getNormalTotalTractions(I);
                ref_fric_I = FricPtr->get_reference_friction_coefficients_(I);
                a_I = FricPtr->get_a_values(I);
                b_I = FricPtr->get_b_values(I);
                ref_slip_rate_I = FricPtr->get_reference_velocities(I);
                average_slope3_eq4[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc3_(ref_slip_rate_I, a_I, b_I, theta_I, dc_I, ref_fric_I, sigma_n_I, press_I, slip_rate_I, ElastMatTimesDDsRates_I, perm_I, fluid_compressibility, viscosity, hydraulic_aperture_I, void_compressibility, fault_porosity, laplacian_pressure[I], Q_i, loading_rate));
                Solution_.setSlipRateAlongSlipDirection(I, previous_DDs_Rates[3 * I + slip_direction] + ((1932. / 2197.) * average_slope1_eq4[I]) - ((7200. / 2197.) * average_slope2_eq4[I]) + ((7296. / 2197.) * average_slope3_eq4[I]), slip_direction);
            }
            for (arma::uword be = 0; be < boundary_elts.n_elem; ++be) {
                apply_pressure_bc(boundary_elts[be]);
            }
            laplacian_pressure = EQSim::CalculatePressureLaplacianUniformMeshOnly(neigh_elts, FaultMesh, Solution_.getPressure(), neumann_bc);
            if (use_hmatrix) {
                // Use the H-matrix multiplication
                // ElastMatTimesDDsRates = hmat.multiply(Solution_.getDDsRates() + Vload_vec);
                ElastMatTimesDDsRates = hmat.multiply(Solution_.getDDsRates());
                ElastMatTimesDDsRates -= constant_stress_loading;
            } else {
                ElastMatTimesDDsRates = ElastMat * (Solution_.getDDsRates() - Vload_vec);
            }
            for (arma::uword I = 0; I < Nelts; ++I) {
                perm_I = PermPtr->getPermeabilities(I);
                dc_I = FricPtr->get_state_evolution_distances_(I);
                Q_i = Q_vector[I];
                hydraulic_aperture_I = previous_fault_hydraulic_aperture[I];
                press_I = Solution_.getPressure(I);
                slip_rate_I = Solution_.getSlipRateAlongSlipDirection(I, slip_direction);
                theta_I = previous_state_variables[I] + ((1932. / 2197.) * average_slope1_eq3[I]) - ((7200. / 2197.) * average_slope2_eq3[I]) + ((7296. / 2197.) * average_slope3_eq3[I]);
                plastic_porosity_I = previous_plastic_fault_porosity[I];
                ElastMatTimesDDsRates_I = ElastMatTimesDDsRates[3 * I + slip_direction];
                average_slope4_eq1[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc1_(perm_I, fluid_compressibility, viscosity, hydraulic_aperture_I, void_compressibility, fault_porosity, laplacian_pressure[I], Q_i));
                Solution_.setPressure(I, previous_pressure[I] + ((439. / 216.) * average_slope1_eq1[I]) - (8. * average_slope2_eq1[I]) + ((3680. / 513.) * average_slope3_eq1[I]) - ((845. / 4104.) * average_slope4_eq1[I]));
                //if (InjectionObj.isHoldingPhase()){
                //    Solution_.setPressure( srcIdx, InjectionObj.getTargetPressure());
                //    laplacian_pressure[srcIdx] = 0.0;}
                average_slope4_eq2[I] = Solution_.getTimeStep() * slip_rate_I;
                average_slope4_eq3[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc2_(theta_I, dc_I, slip_rate_I));
                sigma_n_I = Solution_.getNormalTotalTractions(I);
                ref_fric_I = FricPtr->get_reference_friction_coefficients_(I);
                a_I = FricPtr->get_a_values(I);
                b_I = FricPtr->get_b_values(I);
                ref_slip_rate_I = FricPtr->get_reference_velocities(I);
                average_slope4_eq4[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc3_(ref_slip_rate_I, a_I, b_I, theta_I, dc_I, ref_fric_I, sigma_n_I, press_I, slip_rate_I, ElastMatTimesDDsRates_I, perm_I, fluid_compressibility, viscosity, hydraulic_aperture_I, void_compressibility, fault_porosity, laplacian_pressure[I], Q_i, loading_rate));
                Solution_.setSlipRateAlongSlipDirection(I, previous_DDs_Rates[3 * I + slip_direction] + ((439. / 216.) * average_slope1_eq4[I]) - (8. * average_slope2_eq4[I]) + ((3680. / 513.) * average_slope3_eq4[I]) - ((845. / 4104.) * average_slope4_eq4[I]), slip_direction);
            }
            for (arma::uword be = 0; be < boundary_elts.n_elem; ++be) {
                apply_pressure_bc(boundary_elts[be]);
            }
            laplacian_pressure = EQSim::CalculatePressureLaplacianUniformMeshOnly(neigh_elts, FaultMesh, Solution_.getPressure(), neumann_bc);
            if (use_hmatrix) {
                // ElastMatTimesDDsRates = hmat.multiply(Solution_.getDDsRates() + Vload_vec);
                ElastMatTimesDDsRates = hmat.multiply(Solution_.getDDsRates());
                ElastMatTimesDDsRates -= constant_stress_loading;
            } else {
                ElastMatTimesDDsRates = ElastMat * (Solution_.getDDsRates() - Vload_vec);
            }
            for (arma::uword I = 0; I < Nelts; ++I) {
                perm_I = PermPtr->getPermeabilities(I);
                dc_I = FricPtr->get_state_evolution_distances_(I);
                Q_i = Q_vector[I];
                hydraulic_aperture_I = previous_fault_hydraulic_aperture[I];
                press_I = Solution_.getPressure(I);
                slip_rate_I = Solution_.getSlipRateAlongSlipDirection(I, slip_direction);
                theta_I = previous_state_variables[I] + ((439. / 216.) * average_slope1_eq3[I]) - (8. * average_slope2_eq3[I]) + ((3680. / 513.) * average_slope3_eq3[I]) - ((845. / 4104.) * average_slope4_eq3[I]);
                plastic_porosity_I = previous_plastic_fault_porosity[I];
                ElastMatTimesDDsRates_I = ElastMatTimesDDsRates[3 * I + slip_direction];
                average_slope5_eq1[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc1_(perm_I, fluid_compressibility, viscosity, hydraulic_aperture_I, void_compressibility, fault_porosity, laplacian_pressure[I], Q_i));
                Solution_.setPressure(I, previous_pressure[I] - ((8. / 27.) * average_slope1_eq1[I]) + (2. * average_slope2_eq1[I]) - ((3544. / 2565.) * average_slope3_eq1[I]) + ((1859. / 4104.) * average_slope4_eq1[I]) - ((11. / 40.) * average_slope5_eq1[I]));
                //if (InjectionObj.isHoldingPhase()){
                //    Solution_.setPressure( srcIdx, InjectionObj.getTargetPressure());
                //    laplacian_pressure[srcIdx] = 0.0;}
                average_slope5_eq2[I] = Solution_.getTimeStep() * slip_rate_I;
                average_slope5_eq3[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc2_(theta_I, dc_I, slip_rate_I));
                sigma_n_I = Solution_.getNormalTotalTractions(I);
                ref_fric_I = FricPtr->get_reference_friction_coefficients_(I);
                a_I = FricPtr->get_a_values(I);
                b_I = FricPtr->get_b_values(I);
                ref_slip_rate_I = FricPtr->get_reference_velocities(I);
                average_slope5_eq4[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc3_(ref_slip_rate_I, a_I, b_I, theta_I, dc_I, ref_fric_I, sigma_n_I, press_I, slip_rate_I, ElastMatTimesDDsRates_I, perm_I, fluid_compressibility, viscosity, hydraulic_aperture_I, void_compressibility, fault_porosity, laplacian_pressure[I], Q_i, loading_rate));
                Solution_.setSlipRateAlongSlipDirection(I, previous_DDs_Rates[3 * I + slip_direction] - ((8. / 27.) * average_slope1_eq4[I]) + (2. * average_slope2_eq4[I]) - ((3544. / 2565.) * average_slope3_eq4[I]) + ((1859. / 4104.) * average_slope4_eq4[I]) - ((11. / 40.) * average_slope5_eq4[I]), slip_direction);
            }
            for (arma::uword be = 0; be < boundary_elts.n_elem; ++be) {
                apply_pressure_bc(boundary_elts[be]);
            }
            laplacian_pressure = EQSim::CalculatePressureLaplacianUniformMeshOnly(neigh_elts, FaultMesh, Solution_.getPressure(), neumann_bc);
            if (use_hmatrix) {
                // Use the H-matrix multiplication
                // ElastMatTimesDDsRates = hmat.multiply(Solution_.getDDsRates() + Vload_vec);
                ElastMatTimesDDsRates = hmat.multiply(Solution_.getDDsRates());
                ElastMatTimesDDsRates -= constant_stress_loading;
            } else {
                ElastMatTimesDDsRates = ElastMat * (Solution_.getDDsRates() - Vload_vec);
            }
            for (arma::uword I = 0; I < Nelts; ++I) {
                perm_I = PermPtr->getPermeabilities(I);
                dc_I = FricPtr->get_state_evolution_distances_(I);
                Q_i = Q_vector[I];
                hydraulic_aperture_I = previous_fault_hydraulic_aperture[I];
                press_I = Solution_.getPressure(I);
                slip_rate_I = Solution_.getSlipRateAlongSlipDirection(I, slip_direction);
                theta_I = previous_state_variables[I] - ((8. / 27.) * average_slope1_eq3[I]) + (2. * average_slope2_eq3[I]) - ((3544. / 2565.) * average_slope3_eq3[I]) + ((1859. / 4104.) * average_slope4_eq3[I]) - ((11. / 40.) * average_slope5_eq3[I]);
                plastic_porosity_I = previous_plastic_fault_porosity[I];
                ElastMatTimesDDsRates_I = ElastMatTimesDDsRates[3 * I + slip_direction];
                average_slope6_eq1[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc1_(perm_I, fluid_compressibility, viscosity, hydraulic_aperture_I, void_compressibility, fault_porosity, laplacian_pressure[I], Q_i));
                average_slope6_eq2[I] = Solution_.getTimeStep() * slip_rate_I;
                average_slope6_eq3[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc2_(theta_I, dc_I, slip_rate_I));
                sigma_n_I = Solution_.getNormalTotalTractions(I);
                ref_fric_I = FricPtr->get_reference_friction_coefficients_(I);
                a_I = FricPtr->get_a_values(I);
                b_I = FricPtr->get_b_values(I);
                ref_slip_rate_I = FricPtr->get_reference_velocities(I);
                average_slope6_eq4[I] = Solution_.getTimeStep() * (rhsODEs_->RhsFunc3_(ref_slip_rate_I, a_I, b_I, theta_I, dc_I, ref_fric_I, sigma_n_I, press_I, slip_rate_I, ElastMatTimesDDsRates_I, perm_I, fluid_compressibility, viscosity, hydraulic_aperture_I, void_compressibility, fault_porosity, laplacian_pressure[I], Q_i, loading_rate));
                new_pressure_order4[I] = previous_pressure[I] + ((25. / 216.) * average_slope1_eq1[I]) + ((1408. / 2565.) * average_slope3_eq1[I]) + ((2197. / 4104.) * average_slope4_eq1[I]) - ((1. / 5.) * average_slope5_eq1[I]);
                new_pressure_order5[I] = previous_pressure[I] + ((16. / 135.) * average_slope1_eq1[I]) + ((6656. / 12825.) * average_slope3_eq1[I]) + ((28561. / 56430.) * average_slope4_eq1[I]) - ((9. / 50.) * average_slope5_eq1[I]) + ((2. / 55.) * average_slope6_eq1[I]);
                // Set new slip
                new_slip_order4[I] = previous_DDs[3 * I + slip_direction] + ((25. / 216.) * average_slope1_eq2[I]) + ((1408. / 2565.) * average_slope3_eq2[I]) + ((2197. / 4104.) * average_slope4_eq2[I]) - ((1. / 5.) * average_slope5_eq2[I]);
                new_slip_order5[I] = previous_DDs[3 * I + slip_direction] + ((16. / 135.) * average_slope1_eq2[I]) + ((6656. / 12825.) * average_slope3_eq2[I]) + ((28561. / 56430.) * average_slope4_eq2[I]) - ((9. / 50.) * average_slope5_eq2[I]) + ((2. / 55.) * average_slope6_eq2[I]);
                // Set new state variable
                new_theta_order4[I] = previous_state_variables[I] + ((25. / 216.) * average_slope1_eq3[I]) + ((1408. / 2565.) * average_slope3_eq3[I]) + ((2197. / 4104.) * average_slope4_eq3[I]) - ((1. / 5.) * average_slope5_eq3[I]);
                new_theta_order5[I] = previous_state_variables[I] + ((16. / 135.) * average_slope1_eq3[I]) + ((6656. / 12825.) * average_slope3_eq3[I]) + ((28561. / 56430.) * average_slope4_eq3[I]) - ((9. / 50.) * average_slope5_eq3[I]) + ((2. / 55.) * average_slope6_eq3[I]);
                // Set new plastic fault porosity
                new_plastic_fault_porosity_order4[I] = previous_plastic_fault_porosity[I];
                new_plastic_fault_porosity_order5[I] = previous_plastic_fault_porosity[I];
                // Set new hydraulic aperture
                new_hydraulic_aperture_order4[I] = previous_fault_hydraulic_aperture[I];
                new_hydraulic_aperture_order5[I] = previous_fault_hydraulic_aperture[I];
                // Set new slip rates
                new_slip_rates_order4[I] = previous_DDs_Rates[3 * I + slip_direction] + ((25. / 216.) * average_slope1_eq4[I]) + ((1408. / 2565.) * average_slope3_eq4[I]) + ((2197. / 4104.) * average_slope4_eq4[I]) - ((1. / 5.) * average_slope5_eq4[I]);
                new_slip_rates_order5[I] = previous_DDs_Rates[3 * I + slip_direction] + ((16. / 135.) * average_slope1_eq4[I]) + ((6656. / 12825.) * average_slope3_eq4[I]) + ((28561. / 56430.) * average_slope4_eq4[I]) - ((9. / 50.) * average_slope5_eq4[I]) + ((2. / 55.) * average_slope6_eq4[I]);
            }
            for (arma::uword be = 0; be < boundary_elts.n_elem; ++be) {
                apply_pressure_bc_vec(new_pressure_order4, boundary_elts[be]);
                apply_pressure_bc_vec(new_pressure_order5, boundary_elts[be]);
            }
            // Calculate difference between solution of order 5 and the one of order 4
            for (arma::uword I = 0; I < Nelts; ++I) {
                diff_pressure45[I] = new_pressure_order5[I] - new_pressure_order4[I];
                diff_slip_rates45[I] = new_slip_rates_order5[I] - new_slip_rates_order4[I];
                diff_theta45[I] = new_theta_order5[I] - new_theta_order4[I];
            }
            // Compute relative error
            error_pressure45 = arma::norm(diff_pressure45, 1) / (arma::norm(new_pressure_order5, 1)+1e-12);
            // error_slip_rate45 = arma::norm(diff_slip_rates45, 1) / (arma::norm(new_slip_rates_order5, 1)+1e-12);
            double scale_v = std::max(norm(new_slip_rates_order5,1), 1e-8); // pick a floor
            error_slip_rate45 = norm(diff_slip_rates45,1) / scale_v;

            error_theta45 = arma::norm(diff_theta45, 1) / (arma::norm(new_theta_order5, 1)+1e-12);
            max_error = std::max({error_pressure45, error_slip_rate45, error_theta45});
            std::cout << "Maximum Slip rate at the moment is: " << arma::max(new_slip_rates_order4) << std::endl;
            std::cout << "Maximum State value at the moment is: " << arma::max(new_theta_order4) << std::endl;
            std::cout << "Smallest State value at the moment is: " << arma::min(new_theta_order4) << std::endl;
            std::cout << "Relative error slip rate: " << error_slip_rate45 << std::endl;
            std::cout << "Relative error state variable: " << error_theta45 << std::endl;
            if ((error_pressure45 < solver_parameters.tolerance_RK) && (error_slip_rate45 < solver_parameters.tolerance_RK) && (error_theta45 < solver_parameters.tolerance_RK)) {
                std::cout << "Solution converged after " << iter << " iterations!" << std::endl;
                std::cout << "Relative error pressure: " << error_pressure45 << std::endl;
                std::cout << "Relative error slip rate: " << error_slip_rate45 << std::endl;
                std::cout << "Relative error state variable: " << error_theta45 << std::endl;
                Solution_.setPressure(new_pressure_order5);
                //if (InjectionObj.isHoldingPhase()){
                //    Solution_.setPressure( srcIdx, InjectionObj.getTargetPressure());
                //    laplacian_pressure[srcIdx] = 0.0;}
                Solution_.setSlipAlongSlipDirection(new_slip_order5, slip_direction);  
                Solution_.setSlipRateAlongSlipDirection(new_slip_rates_order5, slip_direction);
                Solution_.setStateVariables(new_theta_order5);
                Solution_.setPlasticFaultPorosity(new_plastic_fault_porosity_order5);
                Solution_.setFaultHydraulicAperture(new_hydraulic_aperture_order5);
                for (arma::uword I = 0; I < Nelts; ++I) {
                    auto cumul_incrm_hydraulic_apert = new_hydraulic_aperture_order5[I] - initial_hydraulic_aperture[I];
                    Solution_.setOpening(cumul_incrm_hydraulic_apert, I);
                }
                double new_shear_stress1, new_shear_stress2, new_normal_stress;
                if (use_hmatrix) {
                    // Use the H-matrix multiplication
                    ElastMatTimesDDs = hmat.multiply(Solution_.getDDs());
                } else {
                    ElastMatTimesDDs = ElastMat * Solution_.getDDs();
                }
                for (arma::uword I = 0; I < Nelts; ++I) {
                    new_shear_stress1 = Solution_.getInsituTractions(3 * I) - ElastMatTimesDDs[3 * I];
                    new_shear_stress2 = Solution_.getInsituTractions(3 * I + 1) - ElastMatTimesDDs[3 * I + 1];
                    new_normal_stress = Solution_.getInsituTractions(3 * I + 2) - ElastMatTimesDDs[3 * I + 2];
                    Solution_.setTotalTractions(I, new_shear_stress1, new_shear_stress2, new_normal_stress);
                }
                Solution_.setCurrentTime(time + Solution_.getTimeStep());
                if ((iter == 1) && ((max_error - solver_parameters.tolerance_RK) > -10e-3)) {
                    /*
                    double dt_ratio = std::numeric_limits<double>::infinity();
                    for (arma::uword i = 0; i < Nelts; ++i) {
                        double sr = std::abs(new_slip_rates_order5[i]);
                        sr = std::max(sr, 1e-12);
                        double dc = FricPtr->get_state_evolution_distances_(i);
                        dt_ratio = std::min(dt_ratio, dc / sr);
                    }
                    double dt_limit = solver_parameters.adaptive_dt_constant * dt_ratio;
                    double amplified_dt = solver_parameters.time_step_amplification_factor * Solution_.getTimeStep();
                    Solution_.setTimeStep(std::min(amplified_dt, dt_limit));*/
                    Solution_.setTimeStep(std::min(solver_parameters.time_step_amplification_factor * Solution_.getTimeStep(), 500.0));
                }
            } else {
                Solution_.setTimeStep(Solution_.getTimeStep() / solver_parameters.time_step_reduction_factor);
                Solution_.setPressure(previous_pressure);
                Solution_.setSlipRateAlongSlipDirection(previous_DDs_Rates, slip_direction);
            }
        }
        //Initiate netcdf file
        if (counter == 1){
            ncfilenameSolution = Solution_.getResPath() + Solution_.getBaseFileName();
            Solution_.initializeNetCDF(ncfilenameSolution);
        }
        //Export the output file
        if ((counter == 1) || (counter % solver_parameters.export_current_solution_every_i_time_steps_ == 0)) {
            json_export_solution = Solution_.createJsonObjectToExport(export_background_stresses);
            filenameSolution = Solution_.getResPath() + Solution_.getBaseFileName() + std::string{"-Solution_at_time_"} + std::to_string(Solution_.getCurrentTime());
            Solution_.exportToJson(json_export_solution, filenameSolution);
            Solution_.exportBaseSolutionToNetCDF(Solution_.getCurrentTime());
        }
        if (counter % solver_parameters.export_stress_every_i_time_steps_ == 0) {
            std::vector<double> slipStdStep = arma::conv_to<std::vector<double>>::from(Solution_.getDDs());
            std::string filenameStress = Solution_.getResPath() + "stress_at_time_" + std::to_string(Solution_.getCurrentTime()) + ".csv";
            stressExporter.ComputeAndExport(slipStdStep, filenameStress);}
        time = Solution_.getCurrentTime();
    }
    //std::vector<double> slipStd = arma::conv_to<std::vector<double>>::from(Solution_.getDDs()); 
    //stressExporter.ComputeAndExport(slipStd,"stress.csv");
    std::vector<double> slipStd = arma::conv_to<std::vector<double>>::from(Solution_.getDDs());
    std::string filenameStressEnd = Solution_.getResPath() + "stress_at_time_" + std::to_string(Solution_.getCurrentTime()) + ".csv";
    stressExporter.ComputeAndExport(slipStd, filenameStressEnd);
    Solution_.closeNetCDF();
}
} // namespace EQSim
