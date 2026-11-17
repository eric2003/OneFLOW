// ==================== src/cfd_solver.cpp (重构版) ====================
#include "cfd_solver.hpp"

namespace cfd {

void CfdSolver::initialize(const CfdConfig& config) {
    config_ = config;
    config_.validate();
    
    // Initialize domain
    domain_ = ComputationalDomain(config_);
    
    // Initialize solution
    solution_ = Solution(domain_);
    
    // Create components using factories
    boundary_condition_ = BoundaryFactory::create_boundary(config_);
    
    // Default: step function for 1D convection
    std::vector<Real> step_params = {1.0, 2.0, 0.5, 1.0};
    initial_condition_ = InitialConditionFactory::create_initial_condition("step", step_params);
    
    // Create residual calculator
    residual_calculator_ = ResidualCalculatorFactory::create_from_config(config_, true); // compact mode
    
    // Create time integrator
    time_integrator_ = TimeIntegratorFactory::create_integrator(
        config_.rk_order,
        std::make_unique<ConvectionResidualCalculator>(
            ReconstructorFactory::create_reconstructor(config_),
            FluxFactory::create_flux_calculator(config_.flux_type),
            BoundaryFactory::create_boundary(config_),
            config_.wave_speed
        )
    );
    
    // Initialize solution
    initialize_solution();
    
    // Adjust time step based on CFL condition
    calculate_dt();
    
    if (config_.verbose) {
        std::cout << "CFD solver initialized:" << std::endl;
        std::cout << "  Time integrator: " << time_integrator_->name() << std::endl;
        std::cout << "  Residual calculator: " << residual_calculator_->name() << std::endl;
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
    // Update residual calculator with new boundary condition
    if (residual_calculator_) {
        residual_calculator_ = ResidualCalculatorFactory::create_convection_residual(
            std::make_unique<EnoReconstructor>(config_.spatial_order),
            FluxFactory::create_flux_calculator(config_.flux_type),
            std::make_unique<PeriodicBoundary>(),
            config_.wave_speed,
            true
        );
    }
    // Re-apply boundary conditions
    boundary_condition_->apply(solution_.u(), domain_);
}

void CfdSolver::set_boundary_condition(const std::string& bc_type, 
                                      Real left_val, Real right_val) {
    set_boundary_condition(BoundaryFactory::create_boundary(bc_type, left_val, right_val));
}

void CfdSolver::set_initial_condition(std::unique_ptr<InitialCondition> ic) {
    initial_condition_ = std::move(ic);
    initialize_solution();
}

void CfdSolver::set_initial_condition(const std::string& ic_type,
                                     const std::vector<Real>& params) {
    set_initial_condition(InitialConditionFactory::create_initial_condition(ic_type, params));
}

void CfdSolver::set_time_integrator(std::unique_ptr<TimeIntegrator> integrator) {
    time_integrator_ = std::move(integrator);
}

void CfdSolver::set_time_integrator(const std::string& method) {
    time_integrator_ = TimeIntegratorFactory::create_integrator(
        method,
        std::move(residual_calculator_)
    );
}

void CfdSolver::set_residual_calculator(std::unique_ptr<ResidualCalculator> calculator) {
    residual_calculator_ = std::move(calculator);
}

void CfdSolver::advance_one_step(Real dt) {
    if (!time_integrator_) {
        throw std::runtime_error("Time integrator not initialized");
    }
    
    // Perform one time step
    time_integrator_->step(solution_.u(), dt, domain_);
    
    // Save old solution
    solution_.copy_to_old();
}

Vector CfdSolver::run_simulation(Real final_time) {
    // Use config final time if not specified
    if (final_time < 0.0) {
        final_time = config_.final_time;
    }
    
    // Ensure solver is properly initialized
    if (!time_integrator_ || !residual_calculator_ || 
        !boundary_condition_ || !initial_condition_) {
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
        std::cout << "  Time integrator: " << time_integrator_->name() << std::endl;
        std::cout << "  Residual calculator: " << residual_calculator_->name() << std::endl;
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
    std::cout << "  Time integrator: " << time_integrator_->name() << " (order " 
              << time_integrator_->order() << ")" << std::endl;
    std::cout << "  Residual calculator: " << residual_calculator_->name() << std::endl;
    std::cout << "  Boundary condition: " << boundary_condition_->name() << std::endl;
    std::cout << "  Initial condition: " << initial_condition_->name() << std::endl;
    std::cout << "==========================================" << std::endl;
}

// Static methods
Vector CfdSolver::run_single_simulation(const CfdConfig& config) {
    CfdSolver solver(config);
    return solver.run_simulation();
}

} // namespace cfd