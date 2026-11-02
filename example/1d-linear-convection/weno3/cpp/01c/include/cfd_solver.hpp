// ==================== include/cfd_solver.hpp ====================
#ifndef CFD_SOLVER_HPP
#define CFD_SOLVER_HPP

#include "config.hpp"
#include "mesh.hpp"
#include "domain.hpp"
#include "solution.hpp"
#include "flux.hpp"
#include "reconstructor.hpp"
#include <memory>
#include <functional>
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
    std::unique_ptr<Reconstructor> reconstructor_;
    std::unique_ptr<FluxCalculator> flux_calculator_;
    
    // Time integration methods
    void runge_kutta_1(Real dt);
    void runge_kutta_2(Real dt);
    
    // Helper methods
    void calculate_dt();
    void compute_residual();
    
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
    
    // Initialize field with custom function
    void initialize_field(const std::function<Real(Real)>& init_func);
    
    // Initialize field with step function (default)
    void initialize_step_function();
    
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
    // Static methods
    // ===================================================================
    static Vector run_single_simulation(const CfdConfig& config);
    static void perform_eno_weno_analysis();
};

} // namespace cfd

#endif // CFD_SOLVER_HPP