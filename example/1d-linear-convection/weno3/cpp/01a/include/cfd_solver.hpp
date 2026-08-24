// ==================== include/cfd_solver.hpp ====================
#ifndef CFD_SOLVER_HPP
#define CFD_SOLVER_HPP

#include "config.hpp"  // 包含配置头文件
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

// Forward declarations
class Mesh;
class ComputationalDomain;
class Reconstructor;
class Solution;

// ===================================================================
// Mesh class
// ===================================================================
class Mesh {
public:
    Real xmin = 0.0;
    Real xmax = 2.0;
    int ncells = 40;
    int nnodes = 0;
    int nx = 0;
    Real L = 0.0;
    Real dx = 0.0;
    Vector x;
    Vector xcc;
    
    Mesh() = default;
    
    // Constructor with CfdConfig
    explicit Mesh(const CfdConfig& config) 
        : xmin(config.xmin), xmax(config.xmax), ncells(config.ncells) {
        init();
    }
    
    // Constructor with explicit parameters
    Mesh(Real xmin_, Real xmax_, int ncells_) 
        : xmin(xmin_), xmax(xmax_), ncells(ncells_) {
        init();
    }
    
    void init() {
        nnodes = ncells + 1;
        nx = ncells;
        L = xmax - xmin;
        dx = L / static_cast<Real>(ncells);
        
        x.resize(nnodes);
        xcc.resize(ncells);
        
        // Node coordinates
        for (int i = 0; i < nnodes; ++i) {
            x[i] = xmin + static_cast<Real>(i) * dx;
        }
        
        // Cell center coordinates
        for (int i = 0; i < ncells; ++i) {
            xcc[i] = 0.5 * (x[i] + x[i + 1]);
        }
        
        if (false) {  // Change to true for debug output
            std::cout << "Mesh initialized:" << std::endl;
            std::cout << "  ncells = " << ncells << std::endl;
            std::cout << "  dx = " << dx << std::endl;
            std::cout << "  L = " << L << std::endl;
        }
    }
    
    void print_info() const {
        std::cout << "Mesh Information:" << std::endl;
        std::cout << "  Domain: [" << xmin << ", " << xmax << "]" << std::endl;
        std::cout << "  Cells: " << ncells << std::endl;
        std::cout << "  Nodes: " << nnodes << std::endl;
        std::cout << "  Cell size (dx): " << dx << std::endl;
        std::cout << "  Domain length: " << L << std::endl;
    }
};

// ===================================================================
// Computational domain class
// ===================================================================
class ComputationalDomain {
public:
    Mesh mesh;
    CfdConfig config;
    int nghosts = 0;
    int ist = 0;
    int ied = 0;
    int ntcells = 0;
    
    ComputationalDomain() = default;
    
    ComputationalDomain(const CfdConfig& cfg) {
        init(cfg);
    }
    
    void init(const CfdConfig& cfg) {
        config = cfg;
        config.validate();
        
        // Initialize mesh from config
        mesh = Mesh(config);
        
        // Calculate number of ghost cells
        if (config.recon_scheme == "eno") {
            nghosts = config.spatial_order;
        } else if (config.recon_scheme == "weno") {
            nghosts = config.spatial_order / 2 + 1;
        } else {
            throw std::runtime_error("Unknown reconstruction scheme: " + config.recon_scheme);
        }
        
        ist = nghosts + 1;
        ied = ist + mesh.ncells - 1;
        ntcells = mesh.ncells + 2 * nghosts;
        
        if (config.verbose) {
            std::cout << "Domain initialized:" << std::endl;
            std::cout << "  mesh.ncells = " << mesh.ncells << std::endl;
            std::cout << "  spatial_order = " << config.spatial_order << std::endl;
            std::cout << "  nghosts = " << nghosts << std::endl;
            std::cout << "  ist = " << ist << ", ied = " << ied << std::endl;
            std::cout << "  ntcells = " << ntcells << std::endl;
        }
    }
    
    void print_info() const {
        std::cout << "Computational Domain:" << std::endl;
        std::cout << "  Total cells (with ghosts): " << ntcells << std::endl;
        std::cout << "  Ghost cells: " << nghosts << " (left and right)" << std::endl;
        std::cout << "  Interior range: [" << ist << ", " << ied << "]" << std::endl;
        std::cout << "  Interior cells: " << (ied - ist + 1) << std::endl;
    }
};

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
public:
    ComputationalDomain domain;
    Vector q_face_left;
    Vector q_face_right;
    Vector flux;
    Vector res;
    Vector u;
    Vector un;
    
    Solution() = default;
    explicit Solution(const ComputationalDomain& dom) { init(dom); }
    
    void init(const ComputationalDomain& dom) {
        domain = dom;
        
        q_face_left.resize(dom.mesh.nnodes, 0.0);
        q_face_right.resize(dom.mesh.nnodes, 0.0);
        flux.resize(dom.mesh.nnodes, 0.0);
        res.resize(dom.mesh.ncells, 0.0);
        u.resize(dom.ntcells, 0.0);
        un.resize(dom.ntcells, 0.0);
    }
    
    void copy_old_field() {
        std::copy(u.begin(), u.end(), un.begin());
    }
    
    // Get interior solution (without ghost cells)
    Vector get_interior_solution() const {
        Vector interior(domain.mesh.ncells);
        std::copy(u.begin() + domain.ist,
                 u.begin() + domain.ist + domain.mesh.ncells,
                 interior.begin());
        return interior;
    }
    
    // Get cell center coordinates
    const Vector& get_cell_centers() const {
        return domain.mesh.xcc;
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
    void apply_periodic_boundary(Vector& u) const;
    void apply_boundary(Vector& u) { apply_periodic_boundary(u); }
    
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