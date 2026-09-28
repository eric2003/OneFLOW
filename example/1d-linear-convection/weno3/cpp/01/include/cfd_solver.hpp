// ==================== include/cfd_solver.hpp ====================
#ifndef CFD_SOLVER_HPP
#define CFD_SOLVER_HPP

#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdexcept>

namespace cfd {

// 类型别名
using Real = double;
using Vector = std::vector<Real>;

// ===================================================================
// 配置结构体
// ===================================================================
struct CfdConfig {
    std::string recon_scheme = "eno";
    int flux_type = 0;
    int rk_order = 1;
    int spatial_order = 3;
    Real wave_speed = 1.0;
    Real final_time = 0.625;
    Real dt = 0.025;
    Real cfl = 0.5;
    
    void print() const {
        std::cout << "CFD Configuration:" << std::endl;
        std::cout << "  recon_scheme: " << recon_scheme << std::endl;
        std::cout << "  flux_type: " << flux_type << std::endl;
        std::cout << "  rk_order: " << rk_order << std::endl;
        std::cout << "  spatial_order: " << spatial_order << std::endl;
        std::cout << "  wave_speed: " << wave_speed << std::endl;
        std::cout << "  final_time: " << final_time << std::endl;
        std::cout << "  dt: " << dt << std::endl;
        std::cout << "  cfl: " << cfl << std::endl;
    }
};

// ===================================================================
// 网格类
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
    Mesh(int n) : ncells(n) { init(); }
    
    void init() {
        nnodes = ncells + 1;
        nx = ncells;
        L = xmax - xmin;
        dx = L / static_cast<Real>(ncells);
        
        x.resize(nnodes);
        xcc.resize(ncells);
        
        // 节点坐标
        for (int i = 0; i < nnodes; ++i) {
            x[i] = xmin + i * dx;
        }
        
        // 单元中心坐标
        for (int i = 0; i < ncells; ++i) {
            xcc[i] = 0.5 * (x[i] + x[i + 1]);
        }
        
        std::cout << "Mesh initialized:" << std::endl;
        std::cout << "  ncells = " << ncells << std::endl;
        std::cout << "  dx = " << dx << std::endl;
        std::cout << "  L = " << L << std::endl;
    }
};

// ===================================================================
// 计算域类
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
    ComputationalDomain(const Mesh& m, const CfdConfig& cfg) {
        init(m, cfg);
    }
    
    void init(const Mesh& m, const CfdConfig& cfg) {
        mesh = m;
        config = cfg;
        
        // 计算鬼单元数
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
        
        std::cout << "Domain initialized:" << std::endl;
        std::cout << "  mesh.ncells = " << mesh.ncells << std::endl;
        std::cout << "  spatial_order = " << config.spatial_order << std::endl;
        std::cout << "  nghosts = " << nghosts << std::endl;
        std::cout << "  ist = " << ist << ", ied = " << ied << std::endl;
        std::cout << "  ntcells = " << ntcells << std::endl;
    }
};

// ===================================================================
// 重构器基类（抽象类）
// ===================================================================
class Reconstructor {
public:
    virtual ~Reconstructor() = default;
    virtual void reconstruct(const Vector& q, Vector& q_face_left, Vector& q_face_right,
                           const ComputationalDomain& domain) = 0;
    virtual std::string name() const = 0;
};

// ===================================================================
// ENO重构器
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
};

// ===================================================================
// WENO重构器
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
};

// ===================================================================
// 解类
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
};

// ===================================================================
// 主CFD求解器类
// ===================================================================
class CfdSolver {
private:
    CfdConfig config_;
    ComputationalDomain domain_;
    Solution solution_;
    std::unique_ptr<Reconstructor> reconstructor_;
    
    // 辅助函数
    void calculate_dt();
    void apply_periodic_boundary(Vector& u) const;
    void apply_boundary(Vector& u) { apply_periodic_boundary(u); }
    
    Real initial_condition(Real x) const;
    Real analytical_solution(Real x, Real t, Real a, Real L) const;
    
    void rusanov_flux(const Vector& qL, const Vector& qR, Vector& flux) const;
    void engquist_osher_flux(const Vector& qL, const Vector& qR, Vector& flux) const;
    void inviscid_flux(const Vector& qL, const Vector& qR, Vector& flux) const;
    
    void residual(const Vector& q);
    
    void runge_kutta_1();
    void runge_kutta_2();
    void runge_kutta();
    
public:
    CfdSolver() = default;
    CfdSolver(const CfdConfig& config, const ComputationalDomain& domain);
    
    void init(const CfdConfig& config, const ComputationalDomain& domain);
    
    void init_field();
    Vector run_simulation(Real final_time);
    
    const Solution& get_solution() const { return solution_; }
    const ComputationalDomain& get_domain() const { return domain_; }
    const CfdConfig& get_config() const { return config_; }
    
    static void perform_eno_weno_analysis();
};

} // namespace cfd

#endif // CFD_SOLVER_HPP