// ==================== include/cfd_solver.hpp ====================
#ifndef CFD_SOLVER_HPP
#define CFD_SOLVER_HPP

#include "config.hpp"
#include "mesh.hpp"
#include "domain.hpp"
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <iomanip>
#include <tuple>
#include <functional>

namespace cfd {

// Type aliases
using Real = double;
using Vector = std::vector<Real>;

// ===================================================================
// Reconstructor base class (abstract)
// ===================================================================
class Reconstructor {
public:
    virtual ~Reconstructor() = default;
    virtual void reconstruct(const Vector& q, Vector& q_face_left, Vector& q_face_right,
                           const ComputationalDomain& domain) = 0;
    virtual std::string name() const = 0;
    virtual int get_order() const = 0;
};

// ===================================================================
// ENO reconstructor
// ===================================================================
class EnoReconstructor : public Reconstructor {
private:
    int spatial_order_;
    int ntcells_;
    std::vector<int> lmc_;
    std::vector<std::vector<Real>> coef_;
    std::vector<std::vector<Real>> dd_;
    
    void init_coef(int order, std::vector<std::vector<Real>>& coef);
    
public:
    EnoReconstructor(int spatial_order, int ntcells);
    ~EnoReconstructor() override = default;
    
    void reconstruct(const Vector& q, Vector& q_face_left, Vector& q_face_right,
                    const ComputationalDomain& domain) override;
    
    std::string name() const override { return "ENO" + std::to_string(spatial_order_); }
    int get_order() const override { return spatial_order_; }
};

// ===================================================================
// WENO reconstructor
// ===================================================================
class WenoReconstructor : public Reconstructor {
private:
    static constexpr Real eps_weno = 1.0e-6;
    
    Real wc3L(Real v1, Real v2, Real v3) const;
    Real wc3R(Real v1, Real v2, Real v3) const;
    
public:
    WenoReconstructor() = default;
    ~WenoReconstructor() override = default;
    
    void reconstruct(const Vector& q, Vector& q_face_left, Vector& q_face_right,
                    const ComputationalDomain& domain) override;
    
    std::string name() const override { return "WENO3"; }
    int get_order() const override { return 3; }
};

// ===================================================================
// Solution class
// ===================================================================
class Solution {
private:
    ComputationalDomain domain_;
    
public:
    Vector q_face_left;
    Vector q_face_right;
    Vector flux;
    Vector res;
    Vector u;   // Current solution with ghost cells
    Vector un;  // Old solution with ghost cells
    
    Solution() = default;
    explicit Solution(const ComputationalDomain& domain) { init(domain); }
    
    void init(const ComputationalDomain& domain) {
        domain_ = domain;
        
        q_face_left.resize(domain_.mesh().nnodes(), 0.0);
        q_face_right.resize(domain_.mesh().nnodes(), 0.0);
        flux.resize(domain_.mesh().nnodes(), 0.0);
        res.resize(domain_.mesh().ncells(), 0.0);
        u = domain_.create_field_with_ghosts<Vector>();
        un = domain_.create_field_with_ghosts<Vector>();
    }
    
    const ComputationalDomain& domain() const { return domain_; }
    
    void copy_old_field() {
        std::copy(u.begin(), u.end(), un.begin());
    }
    
    // Get interior solution (without ghost cells)
    Vector get_interior_solution() const {
        return domain_.extract_interior(u);
    }
    
    // Get cell center coordinates
    const Vector& get_cell_centers() const {
        return domain_.mesh().cell_centers();
    }
    
    // Apply boundary conditions
    void apply_boundary() {
        domain_.apply_periodic_boundary(u);
    }
    
    // Initialize with custom function
    void initialize(const std::function<Real(Real)>& init_func) {
        const auto& mesh = domain_.mesh();
        const auto& xcc = mesh.cell_centers();
        
        for (int i = 0; i < mesh.ncells(); ++i) {
            u[domain_.ist() + i] = init_func(xcc[i]);
        }
        
        apply_boundary();
        copy_old_field();
    }
};

// ===================================================================
// Main CFD solver class
// ===================================================================
class CfdSolver {
private:
    CfdConfig config_;
    ComputationalDomain domain_;
    Solution solution_;
    std::unique_ptr<Reconstructor> reconstructor_;
    
    // Helper functions
    void calculate_dt();
    
    Real initial_condition(Real x) const;
    
    void rusanov_flux(const Vector& qL, const Vector& qR, Vector& flux) const;
    void engquist_osher_flux(const Vector& qL, const Vector& qR, Vector& flux) const;
    void inviscid_flux(const Vector& qL, const Vector& qR, Vector& flux) const;
    
    void residual(const Vector& q);
    
    void runge_kutta_1();
    void runge_kutta_2();
    void runge_kutta();
    
public:
    CfdSolver() = default;
    
    // Constructor with configuration
    explicit CfdSolver(const CfdConfig& config) {
        init(config);
    }
    
    void init(const CfdConfig& config);
    
    // Initialize field with step function
    void init_field();
    
    // Initialize field with custom function
    void init_field(const std::function<Real(Real)>& init_func);
    
    // Run simulation to final time
    Vector run_simulation(Real final_time = -1.0);
    
    // Analytical solution
    Real analytical_solution(Real x, Real t) const;
    
    // Compute analytical solution at current time
    Vector compute_analytical_solution(Real t) const;
    
    // Write solution to file
    void write_solution(const std::string& filename, const Vector& sol) const;
    
    // Getters
    const Solution& get_solution() const { return solution_; }
    const ComputationalDomain& get_domain() const { return domain_; }
    const CfdConfig& get_config() const { return config_; }
    
    // Print solver information
    void print_info() const;
    
    // Static analysis function
    static void perform_eno_weno_analysis();
    
    // Run single simulation
    static Vector run_single_simulation(const CfdConfig& config);
};

} // namespace cfd

#endif // CFD_SOLVER_HPP