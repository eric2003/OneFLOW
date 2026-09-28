// ==================== include/cfd_solver.hpp ====================
#ifndef CFD_SOLVER_HPP
#define CFD_SOLVER_HPP

#include "config.hpp"
#include "mesh.hpp"
#include "domain.hpp"
#include "solution.hpp"
#include "flux.hpp"
#include "reconstructor.hpp"
#include "boundary.hpp"      // 新增
#include "initializer.hpp"   // 新增
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
    std::unique_ptr<BoundaryCondition> boundary_condition_;
    std::unique_ptr<InitialCondition> initial_condition_;
    
    // Time integration methods
    void runge_kutta_1(Real dt);
    void runge_kutta_2(Real dt);
    
    // Helper methods
    void calculate_dt();
    void compute_residual();
    
    // Internal initialization
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
    // Boundary condition methods
    // ===================================================================
    void set_boundary_condition(std::unique_ptr<BoundaryCondition> bc);
    void set_boundary_condition(const std::string& bc_type, 
                               Real left_val = 0.0, Real right_val = 0.0);
    const BoundaryCondition& boundary_condition() const { return *boundary_condition_; }
    
    // ===================================================================
    // Initial condition methods
    // ===================================================================
    void set_initial_condition(std::unique_ptr<InitialCondition> ic);
    void set_initial_condition(const std::string& ic_type,
                              const std::vector<Real>& params = {});
    const InitialCondition& initial_condition() const { return *initial_condition_; }
    
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