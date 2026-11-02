// ==================== include/cfd_solver.hpp ====================
#ifndef CFD_SOLVER_HPP
#define CFD_SOLVER_HPP

#include "config.hpp"
#include "mesh.hpp"
#include "domain.hpp"
#include "solution.hpp"
#include "boundary.hpp"
#include "initializer.hpp"
#include "residual.hpp"
#include "time_integrator.hpp"
#include <memory>
#include <fstream>
#include <filesystem>

namespace cfd {

// ===================================================================
// CFD solver class
// ===================================================================
class CfdSolver {
private:
    CfdConfig config_;
    ComputationalDomain domain_;
    Solution solution_;
    
    // Core components
    std::unique_ptr<TimeIntegrator> time_integrator_;
    std::unique_ptr<ResidualCalculator> residual_calculator_;
    std::unique_ptr<BoundaryCondition> boundary_condition_;
    std::unique_ptr<InitialCondition> initial_condition_;
    
    // Helper methods
    void calculate_dt();
    void initialize_solution();
    
public:
    // ===================================================================
    // Constructors
    // ===================================================================
    CfdSolver() = default;
    
    explicit CfdSolver(const CfdConfig& config) {
        initialize(config);
    }
    
    // ===================================================================
    // Initialization
    // ===================================================================
    void initialize(const CfdConfig& config);
    
    // ===================================================================
    // Component setters
    // ===================================================================
    void set_boundary_condition(std::unique_ptr<BoundaryCondition> bc);
    void set_boundary_condition(const std::string& bc_type, 
                               Real left_val = 0.0, Real right_val = 0.0);
    
    void set_initial_condition(std::unique_ptr<InitialCondition> ic);
    void set_initial_condition(const std::string& ic_type,
                              const std::vector<Real>& params = {});
    
    void set_time_integrator(std::unique_ptr<TimeIntegrator> integrator);
    void set_time_integrator(const std::string& method);
    
    void set_residual_calculator(std::unique_ptr<ResidualCalculator> calculator);
    
    // ===================================================================
    // Component getters
    // ===================================================================
    const BoundaryCondition& boundary_condition() const { return *boundary_condition_; }
    const InitialCondition& initial_condition() const { return *initial_condition_; }
    const TimeIntegrator& time_integrator() const { return *time_integrator_; }
    const ResidualCalculator& residual_calculator() const { return *residual_calculator_; }
    
    // ===================================================================
    // Time integration
    // ===================================================================
    void advance_one_step(Real dt);
    Vector run_simulation(Real final_time = -1.0);
    
    // ===================================================================
    // Analytical solution
    // ===================================================================
    Real analytical_solution(Real x, Real t) const;
    Vector compute_analytical_solution(Real t) const;
    
    // ===================================================================
    // File I/O
    // ===================================================================
    void write_solution(const std::string& filename, const Vector& solution) const;
    void write_current_solution(const std::string& filename) const;
    
    // ===================================================================
    // Getters
    // ===================================================================
    const CfdConfig& config() const { return config_; }
    const ComputationalDomain& domain() const { return domain_; }
    const Solution& solution() const { return solution_; }
    
    // ===================================================================
    // Information
    // ===================================================================
    void print_info() const;
    
    // ===================================================================
    // Static methods (只保留必要的)
    // ===================================================================
    static Vector run_single_simulation(const CfdConfig& config);
};

} // namespace cfd

#endif // CFD_SOLVER_HPP