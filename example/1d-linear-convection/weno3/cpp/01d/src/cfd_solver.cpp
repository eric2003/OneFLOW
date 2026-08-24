// ==================== src/cfd_solver.cpp (修复版) ====================
#include "cfd_solver.hpp"

namespace cfd {

    void CfdSolver::initialize(const CfdConfig& config) {
        config_ = config;
        config_.validate();

        // Initialize domain
        domain_ = ComputationalDomain(config_);

        // Initialize solution
        solution_ = Solution(domain_);

        // Create reconstructor
        reconstructor_ = ReconstructorFactory::create_reconstructor(config_);

        // Create flux calculator
        flux_calculator_ = FluxFactory::create_flux_calculator(config_.flux_type);

        // Create boundary condition (default: periodic for convection)
        boundary_condition_ = BoundaryFactory::create_boundary(config_);

        // Create initial condition (default: step function for 1D convection)
        std::vector<Real> step_params = {1.0, 2.0, 0.5, 1.0};  // low, high, start, end
        initial_condition_ = InitialConditionFactory::create_initial_condition("step", step_params);

        // Initialize solution with default conditions
        initialize_solution();

        // Adjust time step based on CFL condition
        calculate_dt();

        if (config_.verbose) {
            std::cout << "CFD solver initialized:" << std::endl;
            std::cout << "  Reconstructor: " << reconstructor_->name() << std::endl;
            std::cout << "  Flux scheme: " << flux_calculator_->name() << std::endl;
            std::cout << "  Boundary condition: " << boundary_condition_->name() << std::endl;
            std::cout << "  Initial condition: " << initial_condition_->name() << std::endl;
            std::cout << "  Time step: " << config_.dt << std::endl;
        }
    }

    void CfdSolver::initialize_solution() {
        // Apply initial condition
        initial_condition_->initialize_with_ghosts(solution_.u(), domain_);

        // Apply boundary condition
        boundary_condition_->apply(solution_.u(), domain_);

        // Save as old solution
        solution_.copy_to_old();
    }

    void CfdSolver::calculate_dt() {
        // CFL condition: dt <= CFL * dx / |wave_speed|
        Real dt_cfl = config_.cfl * domain_.mesh().dx() / std::abs(config_.wave_speed);

        if (config_.dt > dt_cfl) {
            if (config_.verbose) {
                std::cout << "Adjusting time step for stability:" << std::endl;
                std::cout << "  Original dt = " << config_.dt << std::endl;
                std::cout << "  CFL dt = " << dt_cfl << std::endl;
            }
            config_.dt = dt_cfl;
            if (config_.verbose) {
                std::cout << "  Using dt = " << config_.dt << std::endl;
            }
        }
    }

    void CfdSolver::set_boundary_condition(std::unique_ptr<BoundaryCondition> bc) {
        boundary_condition_ = std::move(bc);
        // Re-apply boundary conditions after changing
        boundary_condition_->apply(solution_.u(), domain_);
    }

    void CfdSolver::set_boundary_condition(const std::string& bc_type, 
        Real left_val, Real right_val) {
        boundary_condition_ = BoundaryFactory::create_boundary(bc_type, left_val, right_val);
        // Re-apply boundary conditions after changing
        boundary_condition_->apply(solution_.u(), domain_);
    }

    void CfdSolver::set_initial_condition(std::unique_ptr<InitialCondition> ic) {
        initial_condition_ = std::move(ic);
        // Reinitialize solution
        initialize_solution();
    }

    void CfdSolver::set_initial_condition(const std::string& ic_type,
        const std::vector<Real>& params) {
        initial_condition_ = InitialConditionFactory::create_initial_condition(ic_type, params);
        // Reinitialize solution
        initialize_solution();
    }

    void CfdSolver::compute_residual() {
        // Apply boundary conditions using the boundary condition object
        boundary_condition_->update_ghost_cells(solution_.u(), domain_);

        // Reconstruct face values
        reconstructor_->reconstruct(solution_.u(), 
            solution_.q_face_left(),
            solution_.q_face_right(),
            domain_);

        // Compute fluxes
        flux_calculator_->compute_flux_array(
            solution_.q_face_left(),
            solution_.q_face_right(),
            solution_.flux(),
            config_.wave_speed
        );

        // Compute residual (flux divergence)
        Real dx_inv = 1.0 / domain_.mesh().dx();
        for (int i = 0; i < domain_.mesh().ncells(); ++i) {
            solution_.res()[i] = -(solution_.flux()[i+1] - solution_.flux()[i]) * dx_inv;
        }
    }

    void CfdSolver::runge_kutta_1(Real dt) {
        // Compute residual
        compute_residual();

        // Update solution
        for (int i = domain_.ist(); i <= domain_.ied(); ++i) {
            int j = i - domain_.ist();
            solution_.u()[i] = solution_.u()[i] + dt * solution_.res()[j];
        }

        // Apply boundary conditions
        boundary_condition_->update_ghost_cells(solution_.u(), domain_);

        // Save old solution
        solution_.copy_to_old();
    }

    void CfdSolver::runge_kutta_2(Real dt) {
        // Stage 1
        compute_residual();

        for (int i = domain_.ist(); i <= domain_.ied(); ++i) {
            int j = i - domain_.ist();
            solution_.u()[i] = solution_.u()[i] + dt * solution_.res()[j];
        }
        boundary_condition_->update_ghost_cells(solution_.u(), domain_);

        // Stage 2
        compute_residual();
        for (int i = domain_.ist(); i <= domain_.ied(); ++i) {
            int j = i - domain_.ist();
            solution_.u()[i] = 0.5 * solution_.un()[i] + 
                0.5 * solution_.u()[i] + 
                0.5 * dt * solution_.res()[j];
        }
        boundary_condition_->update_ghost_cells(solution_.u(), domain_);

        solution_.copy_to_old();
    }

    void CfdSolver::advance_one_step(Real dt) {
        switch (config_.rk_order) {
        case 1:
            runge_kutta_1(dt);
            break;
        case 2:
            runge_kutta_2(dt);
            break;
        default:
            throw std::runtime_error("Unsupported RK order: " + std::to_string(config_.rk_order));
        }
    }

    Vector CfdSolver::run_simulation(Real final_time) {
        // Use config final time if not specified
        if (final_time < 0.0) {
            final_time = config_.final_time;
        }

        // Ensure solution is initialized (safety check)
        if (!boundary_condition_ || !initial_condition_) {
            throw std::runtime_error("Solver not properly initialized. Call initialize() first.");
        }

        Real t = 0.0;
        Real dt_old = config_.dt;
        Real dt = dt_old;
        const int max_steps = 100000;  // Safety limit

        if (config_.verbose) {
            std::cout << "Starting time integration..." << std::endl;
            std::cout << "  Final time: " << final_time << std::endl;
            std::cout << "  Time step: " << dt << std::endl;
            std::cout << "  CFL number: " << config_.cfl << std::endl;
            std::cout << "  Reconstructor: " << reconstructor_->name() << std::endl;
            std::cout << "  Flux scheme: " << flux_calculator_->name() << std::endl;
            std::cout << "  Boundary condition: " << boundary_condition_->name() << std::endl;
            std::cout << "  Initial condition: " << initial_condition_->name() << std::endl;
        }

        int step = 0;
        while (t < final_time - 1.0e-12 && step < max_steps) {
            ++step;

            if (t + dt > final_time) {
                dt = final_time - t;
            }

            advance_one_step(dt);
            t += dt;

            // Progress report
            if (config_.verbose && step % config_.print_interval == 0) {
                std::cout << "  Step " << step << ", Time = " << t << std::endl;
            }
        }

        if (step >= max_steps) {
            std::cout << "Warning: Reached maximum number of steps (" << max_steps << ")" << std::endl;
        }

        config_.dt = dt_old;

        if (config_.verbose) {
            std::cout << "Time integration completed:" << std::endl;
            std::cout << "  Total steps: " << step << std::endl;
            std::cout << "  Final time: " << t << std::endl;
        }

        // Return interior solution (without ghost cells)
        return solution_.get_interior_solution();
    }

    Real CfdSolver::analytical_solution(Real x, Real t) const {
        Real L = domain_.mesh().L();
        Real a = config_.wave_speed;

        // Periodic boundary handling
        Real x_shifted = x - a * t;
        while (x_shifted < 0.0) x_shifted += L;
        while (x_shifted > L) x_shifted -= L;

        // Use the current initial condition to compute analytical solution
        return initial_condition_->evaluate(x_shifted);
    }

    Vector CfdSolver::compute_analytical_solution(Real t) const {
        Vector analytical(domain_.mesh().ncells());
        const auto& xcc = domain_.mesh().cell_centers();

        for (int i = 0; i < domain_.mesh().ncells(); ++i) {
            analytical[i] = analytical_solution(xcc[i], t);
        }
        return analytical;
    }

    void CfdSolver::write_solution(const std::string& filename, const Vector& solution) const {
        solution_.write_to_file(filename, solution, "CFD Solution");
    }

    void CfdSolver::write_current_solution(const std::string& filename) const {
        solution_.write_current_solution(filename, "CFD Current Solution");
    }

    void CfdSolver::print_info() const {
        std::cout << "\nCFD Solver Information:" << std::endl;
        std::cout << "==========================================" << std::endl;
        config_.print_brief();
        std::cout << "\n";
        domain_.print_info();
        std::cout << "  Reconstructor: " << reconstructor_->name() << std::endl;
        std::cout << "  Flux scheme: " << flux_calculator_->name() << std::endl;
        std::cout << "  Boundary condition: " << boundary_condition_->name() << std::endl;
        std::cout << "  Initial condition: " << initial_condition_->name() << std::endl;
        std::cout << "  Time integration: RK" << config_.rk_order << std::endl;
        std::cout << "==========================================" << std::endl;
    }

    // Static methods
    Vector CfdSolver::run_single_simulation(const CfdConfig& config) {
        CfdSolver solver(config);
        // 不再需要调用 initialize_step_function()
        return solver.run_simulation();
    }

    void CfdSolver::perform_eno_weno_analysis() {
        std::cout << "==========================================" << std::endl;
        std::cout << "OneFLOW-CFD Solver: ENO3 vs WENO3 Analysis" << std::endl;
        std::cout << "==========================================" << std::endl;

        // Create configurations
        auto [config_eno, config_weno] = CfdConfig::create_analysis_configs();
        config_eno.output_dir = ".";
        config_weno.output_dir = ".";

        // Run ENO simulation
        std::cout << "\nRunning ENO simulation..." << std::endl;
        config_eno.verbose = true;
        Vector u_eno = run_single_simulation(config_eno);

        // Run WENO simulation
        std::cout << "\nRunning WENO simulation..." << std::endl;
        config_weno.verbose = true;
        Vector u_weno = run_single_simulation(config_weno);

        // Create solver for analytical solution
        CfdSolver solver(config_weno);
        Vector u_analytical = solver.compute_analytical_solution(config_weno.final_time);

        // Write results
        std::cout << "\nWriting results to files..." << std::endl;

        solver.write_solution("eno_results.txt", u_eno);
        solver.write_solution("weno_results.txt", u_weno);
        solver.write_solution("analytical_results.txt", u_analytical);

        // Print statistics
        auto eno_minmax = std::minmax_element(u_eno.begin(), u_eno.end());
        auto weno_minmax = std::minmax_element(u_weno.begin(), u_weno.end());
        auto anal_minmax = std::minmax_element(u_analytical.begin(), u_analytical.end());

        std::cout << "\nSimulation statistics:" << std::endl;
        std::cout << "  ENO min/max: " << *eno_minmax.first << " / " << *eno_minmax.second << std::endl;
        std::cout << "  WENO min/max: " << *weno_minmax.first << " / " << *weno_minmax.second << std::endl;
        std::cout << "  Analytical min/max: " << *anal_minmax.first << " / " << *anal_minmax.second << std::endl;

        std::cout << "\nResults written to:" << std::endl;
        std::cout << "  eno_results.txt" << std::endl;
        std::cout << "  weno_results.txt" << std::endl;
        std::cout << "  analytical_results.txt" << std::endl;

        std::cout << "\nTo generate the comparison plot, run:" << std::endl;
        std::cout << "  python scripts/postprocess.py" << std::endl;
        std::cout << "==========================================" << std::endl;
    }

} // namespace cfd