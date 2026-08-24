// ==================== src/cfd_solver.cpp ====================
#include "cfd_solver.hpp"
#include <algorithm>
#include <iomanip>
#include <sstream>

namespace cfd {

// ===================================================================
// ENO重构器实现
// ===================================================================
EnoReconstructor::EnoReconstructor(int spatial_order, int ntcells)
    : spatial_order_(spatial_order), ntcells_(ntcells) {
    
    lmc_.resize(ntcells, 0);
    
    // 初始化系数矩阵
    coef_.resize(spatial_order + 1);
    for (auto& row : coef_) {
        row.resize(spatial_order, 0.0);
    }
    init_coef(spatial_order, coef_);
    
    // 初始化差商矩阵
    dd_.resize(spatial_order);
    for (auto& row : dd_) {
        row.resize(ntcells, 0.0);
    }
    
    std::cout << "ENO reconstructor initialized:" << std::endl;
    std::cout << "  spatial_order = " << spatial_order_ << std::endl;
    std::cout << "  ntcells = " << ntcells_ << std::endl;
}

void EnoReconstructor::init_coef(int order, std::vector<std::vector<Real>>& coef) {
    switch (order) {
        case 1:
            coef[0][0] = 1.0;
            coef[1][0] = 1.0;
            break;
            
        case 2:
            coef[0][0] = 3.0/2.0;  coef[0][1] = -1.0/2.0;
            coef[1][0] = 1.0/2.0;  coef[1][1] =  1.0/2.0;
            coef[2][0] = -1.0/2.0; coef[2][1] =  3.0/2.0;
            break;
            
        case 3:
            coef[0][0] = 11.0/6.0; coef[0][1] = -7.0/6.0; coef[0][2] =  1.0/3.0;
            coef[1][0] =  1.0/3.0; coef[1][1] =  5.0/6.0; coef[1][2] = -1.0/6.0;
            coef[2][0] = -1.0/6.0; coef[2][1] =  5.0/6.0; coef[2][2] =  1.0/3.0;
            coef[3][0] =  1.0/3.0; coef[3][1] = -7.0/6.0; coef[3][2] = 11.0/6.0;
            break;
            
        default:
            throw std::runtime_error("Unsupported spatial order: " + std::to_string(order));
    }
}

void EnoReconstructor::reconstruct(const Vector& q, Vector& q_face_left,
                                  Vector& q_face_right, const ComputationalDomain& domain) {
    const int nghosts = domain.nghosts;
    const int ist = domain.ist;
    const int ied = domain.ied;
    const int ntcells = domain.ntcells;
    
    // 1. 差商计算
    for (int i = 0; i < ntcells; ++i) {
        dd_[0][i] = q[i];
    }
    
    for (int m = 1; m < spatial_order_; ++m) {
        for (int j = 0; j < ntcells - m; ++j) {
            dd_[m][j] = dd_[m-1][j+1] - dd_[m-1][j];
        }
    }
    
    // 2. 选择最光滑模板
    for (int i = ist - 1; i <= ied; ++i) {
        lmc_[i] = i;
        for (int m = 1; m < spatial_order_; ++m) {
            if (std::abs(dd_[m][lmc_[i] - 1]) < std::abs(dd_[m][lmc_[i]])) {
                lmc_[i] = lmc_[i] - 1;
            }
        }
    }
    
    // 3. 重构界面值
    for (int i = ist; i <= ied; ++i) {
        int j = i - ist;  // 0-based index for interfaces
        
        int k1 = lmc_[i - 1];
        int k2 = lmc_[i];
        int r1 = (i - 1) - k1;
        int r2 = i - k2;
        
        q_face_left[j] = 0.0;
        q_face_right[j] = 0.0;
        
        for (int m = 0; m < spatial_order_; ++m) {
            q_face_left[j] += q[k1 + m] * coef_[r1 + 1][m];
            q_face_right[j] += q[k2 + m] * coef_[r2][m];
        }
    }
}

// ===================================================================
// WENO重构器实现
// ===================================================================
Real WenoReconstructor::wc3L(Real v1, Real v2, Real v3) const {
    // 光滑度指示器
    Real s0 = (v3 - v2) * (v3 - v2);
    Real s1 = (v2 - v1) * (v2 - v1);
    
    // 非线性权值
    Real d0 = 2.0/3.0;
    Real d1 = 1.0/3.0;
    
    Real c0 = d0 / ((eps_weno + s0) * (eps_weno + s0));
    Real c1 = d1 / ((eps_weno + s1) * (eps_weno + s1));
    
    Real w0 = c0 / (c0 + c1);
    Real w1 = c1 / (c0 + c1);
    
    // 候选模板
    Real q0 = 0.5 * v2 + 0.5 * v3;
    Real q1 = -0.5 * v1 + 1.5 * v2;
    
    // 重构值
    return w0 * q0 + w1 * q1;
}

Real WenoReconstructor::wc3R(Real v1, Real v2, Real v3) const {
    // 光滑度指示器
    Real s0 = (v2 - v1) * (v2 - v1);
    Real s1 = (v3 - v2) * (v3 - v2);
    
    // 非线性权值
    Real d0 = 2.0/3.0;
    Real d1 = 1.0/3.0;
    
    Real c0 = d0 / ((eps_weno + s0) * (eps_weno + s0));
    Real c1 = d1 / ((eps_weno + s1) * (eps_weno + s1));
    
    Real w0 = c0 / (c0 + c1);
    Real w1 = c1 / (c0 + c1);
    
    // 候选模板
    Real q0 = 0.5 * v1 + 0.5 * v2;
    Real q1 = 1.5 * v2 - 0.5 * v3;
    
    // 重构值
    return w0 * q0 + w1 * q1;
}

void WenoReconstructor::reconstruct(const Vector& q, Vector& q_face_left,
                                   Vector& q_face_right, const ComputationalDomain& domain) {
    const int nghosts = domain.nghosts;
    const int ist = domain.ist;
    const int ied = domain.ied;
    
    // 左界面重构
    for (int i = ist - 1; i < ied; ++i) {
        int j = i - (ist - 1);
        q_face_left[j] = wc3L(q[i-1], q[i], q[i+1]);
    }
    
    // 右界面重构
    for (int i = ist; i <= ied; ++i) {
        int j = i - ist;
        q_face_right[j] = wc3R(q[i-1], q[i], q[i+1]);
    }
}

// ===================================================================
// CFD求解器实现
// ===================================================================
CfdSolver::CfdSolver(const CfdConfig& config, const ComputationalDomain& domain) {
    init(config, domain);
}

void CfdSolver::init(const CfdConfig& config, const ComputationalDomain& domain) {
    config_ = config;
    domain_ = domain;
    solution_.init(domain);
    
    // 创建重构器
    if (config.recon_scheme == "eno") {
        reconstructor_ = std::make_unique<EnoReconstructor>(config.spatial_order, domain.ntcells);
    } else if (config.recon_scheme == "weno") {
        reconstructor_ = std::make_unique<WenoReconstructor>();
    } else {
        throw std::runtime_error("Unknown reconstruction scheme: " + config.recon_scheme);
    }
    
    // 根据CFL条件调整时间步长
    calculate_dt();
    
    std::cout << "CFD solver initialized with " << reconstructor_->name() << " scheme" << std::endl;
}

void CfdSolver::calculate_dt() {
    // CFL条件: dt <= CFL * dx / |wave_speed|
    Real dt_cfl = config_.cfl * domain_.mesh.dx / std::abs(config_.wave_speed);
    
    if (config_.dt > dt_cfl) {
        std::cout << "Adjusting time step for stability:" << std::endl;
        std::cout << "  Original dt = " << config_.dt << std::endl;
        std::cout << "  CFL dt = " << dt_cfl << std::endl;
        config_.dt = dt_cfl;
        std::cout << "  Using dt = " << config_.dt << std::endl;
    }
}

void CfdSolver::apply_periodic_boundary(Vector& u) const {
    // 复制内部单元到鬼单元
    // 左鬼单元 = 右内部单元
    for (int i = 0; i < domain_.nghosts; ++i) {
        int j = domain_.ist - 1 - i;
        u[j] = u[domain_.ied - domain_.nghosts + i];
    }
    
    // 右鬼单元 = 左内部单元
    for (int i = 0; i < domain_.nghosts; ++i) {
        int j = domain_.ied + i;
        u[j] = u[domain_.ist + i];
    }
}

Real CfdSolver::initial_condition(Real x) const {
    if (0.5 <= x && x <= 1.0) {
        return 2.0;
    } else {
        return 1.0;
    }
}

Real CfdSolver::analytical_solution(Real x, Real t, Real a, Real L) const {
    Real x_shifted = std::fmod(x - a * t + L, L);
    return initial_condition(x_shifted);
}

void CfdSolver::init_field() {
    for (int i = 0; i < domain_.mesh.ncells; ++i) {
        if (0.5 <= domain_.mesh.xcc[i] && domain_.mesh.xcc[i] <= 1.0) {
            solution_.u[domain_.ist + i] = 2.0;
        } else {
            solution_.u[domain_.ist + i] = 1.0;
        }
    }
    
    apply_boundary(solution_.u);
    solution_.copy_old_field();
}

void CfdSolver::rusanov_flux(const Vector& qL, const Vector& qR, Vector& flux) const {
    Real c_L = config_.wave_speed;
    Real c_R = config_.wave_speed;
    
    for (size_t i = 0; i < flux.size(); ++i) {
        Real u_L = qL[i];
        Real u_R = qR[i];
        Real F_L = c_L * u_L;
        Real F_R = c_R * u_R;
        Real Smax = std::max(std::abs(c_L), std::abs(c_R));
        flux[i] = 0.5 * (F_L + F_R) - 0.5 * Smax * (u_R - u_L);
    }
}

void CfdSolver::engquist_osher_flux(const Vector& qL, const Vector& qR, Vector& flux) const {
    Real c = config_.wave_speed;
    Real cp = 0.5 * (c + std::abs(c));
    Real cm = 0.5 * (c - std::abs(c));
    
    for (size_t i = 0; i < flux.size(); ++i) {
        Real u_L = qL[i];
        Real u_R = qR[i];
        flux[i] = cp * u_L + cm * u_R;
    }
}

void CfdSolver::inviscid_flux(const Vector& qL, const Vector& qR, Vector& flux) const {
    if (config_.flux_type == 0) {
        rusanov_flux(qL, qR, flux);
    } else {
        engquist_osher_flux(qL, qR, flux);
    }
}

void CfdSolver::residual(const Vector& q) {
    // 应用边界条件
    Vector q_temp = q;
    apply_boundary(q_temp);
    
    // 重构
    reconstructor_->reconstruct(q_temp, solution_.q_face_left, 
                               solution_.q_face_right, domain_);
    
    // 计算通量
    inviscid_flux(solution_.q_face_left, solution_.q_face_right, solution_.flux);
    
    // 计算残差
    Real dx_inv = 1.0 / domain_.mesh.dx;
    for (int i = 0; i < domain_.mesh.ncells; ++i) {
        solution_.res[i] = -(solution_.flux[i+1] - solution_.flux[i]) * dx_inv;
    }
}

void CfdSolver::runge_kutta_1() {
    Real dt = config_.dt;
    
    // 计算残差
    residual(solution_.u);
    
    // 更新解
    for (int i = domain_.ist; i <= domain_.ied; ++i) {
        int j = i - domain_.ist;
        solution_.u[i] = solution_.u[i] + dt * solution_.res[j];
    }
    
    // 应用边界条件
    apply_boundary(solution_.u);
    
    // 保存旧解
    solution_.copy_old_field();
}

void CfdSolver::runge_kutta_2() {
    Real dt = config_.dt;
    
    // 阶段1
    residual(solution_.u);
    
    for (int i = domain_.ist; i <= domain_.ied; ++i) {
        int j = i - domain_.ist;
        solution_.u[i] = solution_.u[i] + dt * solution_.res[j];
    }
    apply_boundary(solution_.u);
    
    // 阶段2
    residual(solution_.u);
    for (int i = domain_.ist; i <= domain_.ied; ++i) {
        int j = i - domain_.ist;
        solution_.u[i] = 0.5 * solution_.un[i] + 
                         0.5 * solution_.u[i] + 
                         0.5 * dt * solution_.res[j];
    }
    apply_boundary(solution_.u);
    
    solution_.copy_old_field();
}

void CfdSolver::runge_kutta() {
    switch (config_.rk_order) {
        case 1:
            runge_kutta_1();
            break;
        case 2:
            runge_kutta_2();
            break;
        default:
            runge_kutta_1();
            break;
    }
}

Vector CfdSolver::run_simulation(Real final_time) {
    Vector result(domain_.mesh.ncells);
    
    Real t = 0.0;
    Real dt_old = config_.dt;
    Real dt = dt_old;
    const int max_steps = 10000;  // 安全限制
    
    std::cout << "Starting time integration..." << std::endl;
    std::cout << "  Final time: " << final_time << std::endl;
    std::cout << "  Time step: " << dt << std::endl;
    std::cout << "  CFL number: " << config_.cfl << std::endl;
    
    int step = 0;
    while (t < final_time - 1.0e-12 && step < max_steps) {
        ++step;
        
        if (t + dt > final_time) {
            dt = final_time - t;
        }
        
        config_.dt = dt;
        runge_kutta();
        t += dt;
        
        // 进度报告
        if (step % 100 == 0) {
            std::cout << "  Step " << step << ", Time = " << t << std::endl;
        }
    }
    
    if (step >= max_steps) {
        std::cout << "Warning: Reached maximum number of steps (" << max_steps << ")" << std::endl;
    }
    
    config_.dt = dt_old;
    
    std::cout << "Time integration completed:" << std::endl;
    std::cout << "  Total steps: " << step << std::endl;
    std::cout << "  Final time: " << t << std::endl;
    
    // 提取物理解（无鬼单元）
    std::copy(solution_.u.begin() + domain_.ist,
             solution_.u.begin() + domain_.ist + domain_.mesh.ncells,
             result.begin());
    
    return result;
}

void CfdSolver::perform_eno_weno_analysis() {
    std::cout << "==========================================" << std::endl;
    std::cout << "OneFLOW-CFD Solver: ENO3 vs WENO3 Analysis" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // 初始化网格
    Mesh mesh(40);
    mesh.init();
    
    std::cout << "Mesh parameters:" << std::endl;
    std::cout << "  ncells = " << mesh.ncells << std::endl;
    std::cout << "  dx = " << mesh.dx << std::endl;
    std::cout << "  L = " << mesh.L << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // 配置ENO3
    CfdConfig config_eno3;
    config_eno3.recon_scheme = "eno";
    config_eno3.spatial_order = 3;
    config_eno3.flux_type = 0;
    config_eno3.rk_order = 2;
    config_eno3.wave_speed = 1.0;
    config_eno3.final_time = 0.625;
    config_eno3.cfl = 1.0;
    config_eno3.dt = 0.0025;
    
    // 配置WENO3
    CfdConfig config_weno3;
    config_weno3.recon_scheme = "weno";
    config_weno3.spatial_order = 3;
    config_weno3.flux_type = 0;
    config_weno3.rk_order = 2;
    config_weno3.wave_speed = 1.0;
    config_weno3.final_time = 0.625;
    config_weno3.cfl = 1.0;
    config_weno3.dt = 0.0025;
    
    // 创建计算域
    ComputationalDomain domain_eno3(mesh, config_eno3);
    ComputationalDomain domain_weno3(mesh, config_weno3);
    
    // 创建CFD求解器
    CfdSolver solver_eno3(config_eno3, domain_eno3);
    CfdSolver solver_weno3(config_weno3, domain_weno3);
    
    // 运行ENO模拟
    std::cout << "==========================================" << std::endl;
    std::cout << "Running ENO simulation..." << std::endl;
    std::cout << "  Scheme: ENO" << config_eno3.spatial_order << std::endl;
    std::cout << "  Time step: " << config_eno3.dt << std::endl;
    std::cout << "==========================================" << std::endl;
    
    solver_eno3.init_field();
    Vector u_eno = solver_eno3.run_simulation(config_eno3.final_time);
    
    // 运行WENO模拟
    std::cout << "==========================================" << std::endl;
    std::cout << "Running WENO simulation..." << std::endl;
    std::cout << "  Scheme: WENO" << config_weno3.spatial_order << std::endl;
    std::cout << "  Time step: " << config_weno3.dt << std::endl;
    std::cout << "==========================================" << std::endl;
    
    solver_weno3.init_field();
    Vector u_weno = solver_weno3.run_simulation(config_weno3.final_time);
    
    // 计算解析解
    std::cout << "Computing analytical solution..." << std::endl;
    Vector u_analytical(mesh.ncells);
    for (int i = 0; i < mesh.ncells; ++i) {
        u_analytical[i] = solver_weno3.analytical_solution(
            mesh.xcc[i], config_weno3.final_time,
            config_weno3.wave_speed, mesh.L);
    }
    
    // 写入结果文件
    std::cout << "Writing results to files..." << std::endl;
    
    // 写入ENO结果
    std::ofstream eno_file("eno_results.txt");
    eno_file << "# x, u (ENO)\n";
    eno_file << std::fixed << std::setprecision(6);
    for (int i = 0; i < mesh.ncells; ++i) {
        eno_file << std::setw(12) << mesh.xcc[i] 
                << std::setw(12) << u_eno[i] << "\n";
    }
    eno_file.close();
    
    // 写入WENO结果
    std::ofstream weno_file("weno_results.txt");
    weno_file << "# x, u (WENO)\n";
    weno_file << std::fixed << std::setprecision(6);
    for (int i = 0; i < mesh.ncells; ++i) {
        weno_file << std::setw(12) << mesh.xcc[i] 
                 << std::setw(12) << u_weno[i] << "\n";
    }
    weno_file.close();
    
    // 写入解析结果
    std::ofstream analytical_file("analytical_results.txt");
    analytical_file << "# x, u (Analytical)\n";
    analytical_file << std::fixed << std::setprecision(6);
    for (int i = 0; i < mesh.ncells; ++i) {
        analytical_file << std::setw(12) << mesh.xcc[i] 
                       << std::setw(12) << u_analytical[i] << "\n";
    }
    analytical_file.close();
    
    // 打印统计信息
    std::cout << "==========================================" << std::endl;
    std::cout << "Simulation statistics:" << std::endl;
    
    auto [eno_min, eno_max] = std::minmax_element(u_eno.begin(), u_eno.end());
    auto [weno_min, weno_max] = std::minmax_element(u_weno.begin(), u_weno.end());
    auto [anal_min, anal_max] = std::minmax_element(u_analytical.begin(), u_analytical.end());
    
    std::cout << "  ENO min/max: " << *eno_min << " / " << *eno_max << std::endl;
    std::cout << "  WENO min/max: " << *weno_min << " / " << *weno_max << std::endl;
    std::cout << "  Analytical min/max: " << *anal_min << " / " << *anal_max << std::endl;
    
    std::cout << "==========================================" << std::endl;
    std::cout << "Simulation completed successfully!" << std::endl;
    std::cout << "Results written to:" << std::endl;
    std::cout << "  eno_results.txt" << std::endl;
    std::cout << "  weno_results.txt" << std::endl;
    std::cout << "  analytical_results.txt" << std::endl;
    std::cout << std::endl;
    std::cout << "To generate the comparison plot, run:" << std::endl;
    std::cout << "  python scripts/postprocess.py" << std::endl;
    std::cout << "==========================================" << std::endl;
}

} // namespace cfd