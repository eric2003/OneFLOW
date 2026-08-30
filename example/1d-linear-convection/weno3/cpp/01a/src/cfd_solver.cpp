// ==================== src/cfd_solver.cpp ====================
#include "cfd_solver.hpp"
#include <functional>
#include <filesystem>

namespace cfd {

// ===================================================================
// ENO reconstructor implementation
// ===================================================================
EnoReconstructor::EnoReconstructor(int spatial_order, int ntcells)
    : spatial_order_(spatial_order), ntcells_(ntcells) {
    
    if (spatial_order_ < 1 || spatial_order_ > 3) {
        throw std::runtime_error("ENO reconstructor only supports order 1, 2, or 3");
    }
    
    lmc_.resize(ntcells, 0);
    
    // Initialize coefficient matrix
    coef_.resize(spatial_order + 1);
    for (auto& row : coef_) {
        row.resize(spatial_order, 0.0);
    }
    init_coef(spatial_order, coef_);
    
    // Initialize divided difference matrix
    dd_.resize(spatial_order);
    for (auto& row : dd_) {
        row.resize(ntcells, 0.0);
    }
}

void EnoReconstructor::init_coef(int order, std::vector<std::vector<Real>>& coef) {
    // Clear coefficients
    for (auto& row : coef) {
        std::fill(row.begin(), row.end(), 0.0);
    }
    
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
    
    // 1. Compute divided differences
    for (int i = 0; i < ntcells; ++i) {
        dd_[0][i] = q[i];
    }
    
    for (int m = 1; m < spatial_order_; ++m) {
        for (int j = 0; j < ntcells - m; ++j) {
            dd_[m][j] = dd_[m-1][j+1] - dd_[m-1][j];
        }
    }
    
    // 2. Select smoothest stencil
    for (int i = ist - 1; i <= ied; ++i) {
        lmc_[i] = i;
        for (int m = 1; m < spatial_order_; ++m) {
            if (std::abs(dd_[m][lmc_[i] - 1]) < std::abs(dd_[m][lmc_[i]])) {
                lmc_[i] = lmc_[i] - 1;
            }
        }
    }
    
    // 3. Reconstruct face values
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
// WENO reconstructor implementation
// ===================================================================
Real WenoReconstructor::wc3L(Real v1, Real v2, Real v3) const {
    // Smoothness indicators
    Real s0 = (v3 - v2) * (v3 - v2);
    Real s1 = (v2 - v1) * (v2 - v1);
    
    // Nonlinear weights
    Real d0 = 2.0/3.0;
    Real d1 = 1.0/3.0;
    
    Real c0 = d0 / ((eps_weno + s0) * (eps_weno + s0));
    Real c1 = d1 / ((eps_weno + s1) * (eps_weno + s1));
    
    Real w0 = c0 / (c0 + c1);
    Real w1 = c1 / (c0 + c1);
    
    // Candidate stencils
    Real q0 = 0.5 * v2 + 0.5 * v3;
    Real q1 = -0.5 * v1 + 1.5 * v2;
    
    // Reconstructed value
    return w0 * q0 + w1 * q1;
}

Real WenoReconstructor::wc3R(Real v1, Real v2, Real v3) const {
    // Smoothness indicators
    Real s0 = (v2 - v1) * (v2 - v1);
    Real s1 = (v3 - v2) * (v3 - v2);
    
    // Nonlinear weights
    Real d0 = 2.0/3.0;
    Real d1 = 1.0/3.0;
    
    Real c0 = d0 / ((eps_weno + s0) * (eps_weno + s0));
    Real c1 = d1 / ((eps_weno + s1) * (eps_weno + s1));
    
    Real w0 = c0 / (c0 + c1);
    Real w1 = c1 / (c0 + c1);
    
    // Candidate stencils
    Real q0 = 0.5 * v1 + 0.5 * v2;
    Real q1 = 1.5 * v2 - 0.5 * v3;
    
    // Reconstructed value
    return w0 * q0 + w1 * q1;
}

void WenoReconstructor::reconstruct(const Vector& q, Vector& q_face_left,
                                   Vector& q_face_right, const ComputationalDomain& domain) {
    const int ist = domain.ist;
    const int ied = domain.ied;
    
    // Left interface reconstruction
    for (int i = ist - 1; i < ied; ++i) {
        int j = i - (ist - 1);
        q_face_left[j] = wc3L(q[i-1], q[i], q[i+1]);
    }
    
    // Right interface reconstruction
    for (int i = ist; i <= ied; ++i) {
        int j = i - ist;
        q_face_right[j] = wc3R(q[i-1], q[i], q[i+1]);
    }
}

// ===================================================================
// CFD solver implementation
// ===================================================================
void CfdSolver::init(const CfdConfig& config) {
    config_ = config;
    config_.validate();
    
    // Initialize domain
    domain_.init(config_);
    
    // Initialize solution
    solution_.init(domain_);
    
    // Create reconstructor
    if (config_.recon_scheme == "eno") {
        reconstructor_ = std::make_unique<EnoReconstructor>(config_.spatial_order, domain_.ntcells);
    } else if (config_.recon_scheme == "weno") {
        reconstructor_ = std::make_unique<WenoReconstructor>();
    } else {
        throw std::runtime_error("Unknown reconstruction scheme: " + config_.recon_scheme);
    }
    
    // Adjust time step based on CFL condition
    calculate_dt();
    
    if (config_.verbose) {
        std::cout << "CFD solver initialized with " << reconstructor_->name() << " scheme" << std::endl;
    }
}

void CfdSolver::calculate_dt() {
    // CFL condition: dt <= CFL * dx / |wave_speed|
    Real dt_cfl = config_.cfl * domain_.mesh.dx / std::abs(config_.wave_speed);
    
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

void CfdSolver::apply_periodic_boundary(Vector& u) const {
    // Copy interior cells to ghost cells
    // Left ghost cells = right interior cells
    for (int i = 0; i < domain_.nghosts; ++i) {
        int j = domain_.ist - 1 - i;
        u[j] = u[domain_.ied - domain_.nghosts + i];
    }
    
    // Right ghost cells = left interior cells
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

Real CfdSolver::analytical_solution(Real x, Real t) const {
    Real L = domain_.mesh.L;
    Real a = config_.wave_speed;
    Real x_shifted = std::fmod(x - a * t + L, L);
    if (x_shifted < 0.0) x_shifted += L;  // Ensure positive
    return initial_condition(x_shifted);
}

Vector CfdSolver::compute_analytical_solution(Real t) const {
    Vector analytical(domain_.mesh.ncells);
    for (int i = 0; i < domain_.mesh.ncells; ++i) {
        analytical[i] = analytical_solution(domain_.mesh.xcc[i], t);
    }
    return analytical;
}

void CfdSolver::init_field() {
    // Initialize with step function
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

void CfdSolver::init_field(const std::function<Real(Real)>& init_func) {
    // Initialize with custom function
    for (int i = 0; i < domain_.mesh.ncells; ++i) {
        solution_.u[domain_.ist + i] = init_func(domain_.mesh.xcc[i]);
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
    // Apply boundary conditions
    Vector q_temp = q;
    apply_boundary(q_temp);
    
    // Reconstruction
    reconstructor_->reconstruct(q_temp, solution_.q_face_left, 
                               solution_.q_face_right, domain_);
    
    // Compute fluxes
    inviscid_flux(solution_.q_face_left, solution_.q_face_right, solution_.flux);
    
    // Compute residual
    Real dx_inv = 1.0 / domain_.mesh.dx;
    for (int i = 0; i < domain_.mesh.ncells; ++i) {
        solution_.res[i] = -(solution_.flux[i+1] - solution_.flux[i]) * dx_inv;
    }
}

void CfdSolver::runge_kutta_1() {
    Real dt = config_.dt;
    
    // Compute residual
    residual(solution_.u);
    
    // Update solution
    for (int i = domain_.ist; i <= domain_.ied; ++i) {
        int j = i - domain_.ist;
        solution_.u[i] = solution_.u[i] + dt * solution_.res[j];
    }
    
    // Apply boundary conditions
    apply_boundary(solution_.u);
    
    // Save old solution
    solution_.copy_old_field();
}

void CfdSolver::runge_kutta_2() {
    Real dt = config_.dt;
    
    // Stage 1
    residual(solution_.u);
    
    for (int i = domain_.ist; i <= domain_.ied; ++i) {
        int j = i - domain_.ist;
        solution_.u[i] = solution_.u[i] + dt * solution_.res[j];
    }
    apply_boundary(solution_.u);
    
    // Stage 2
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
            throw std::runtime_error("Unsupported RK order: " + std::to_string(config_.rk_order));
    }
}

Vector CfdSolver::run_simulation(Real final_time) {
    // Use config final time if not specified
    if (final_time < 0.0) {
        final_time = config_.final_time;
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
    }
    
    int step = 0;
    while (t < final_time - 1.0e-12 && step < max_steps) {
        ++step;
        
        if (t + dt > final_time) {
            dt = final_time - t;
        }
        
        config_.dt = dt;
        runge_kutta();
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

void CfdSolver::write_solution(const std::string& filename, const Vector& sol) const {
    // Create output directory if it doesn't exist
    std::filesystem::path output_path(config_.output_dir);
    if (!std::filesystem::exists(output_path)) {
        std::filesystem::create_directories(output_path);
    }
    
    std::ofstream file(output_path / filename);
    file << "# x, u\n";
    file << std::fixed << std::setprecision(6);
    
    for (int i = 0; i < domain_.mesh.ncells; ++i) {
        file << std::setw(12) << domain_.mesh.xcc[i] 
             << std::setw(12) << sol[i] << "\n";
    }
    
    file.close();
}

void CfdSolver::print_info() const {
    std::cout << "\nCFD Solver Information:" << std::endl;
    std::cout << "==========================================" << std::endl;
    config_.print_brief();
    std::cout << "\n";
    domain_.print_info();
    std::cout << "  Reconstructor: " << reconstructor_->name() << std::endl;
    std::cout << "  Flux scheme: " << (config_.flux_type == 0 ? "Rusanov" : "Engquist-Osher") << std::endl;
    std::cout << "  Time integration: RK" << config_.rk_order << std::endl;
    std::cout << "==========================================" << std::endl;
}

// Static methods
Vector CfdSolver::run_single_simulation(const CfdConfig& config) {
    CfdSolver solver(config);
    solver.init_field();
    return solver.run_simulation();
}

void CfdSolver::perform_eno_weno_analysis() {
    std::cout << "==========================================" << std::endl;
    std::cout << "OneFLOW-CFD Solver: ENO3 vs WENO3 Analysis" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // Create configurations
    auto [config_eno, config_weno] = CfdConfig::create_analysis_configs();
    
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
    
    // Ensure output directory exists
    std::filesystem::path output_path(config_weno.output_dir);
    if (!std::filesystem::exists(output_path)) {
        std::filesystem::create_directories(output_path);
    }
    
    // Write ENO results
    std::ofstream eno_file(output_path / config_eno.eno_filename);
    eno_file << "# x, u (ENO)\n";
    eno_file << std::fixed << std::setprecision(6);
    for (int i = 0; i < solver.get_domain().mesh.ncells; ++i) {
        eno_file << std::setw(12) << solver.get_domain().mesh.xcc[i] 
                << std::setw(12) << u_eno[i] << "\n";
    }
    eno_file.close();
    
    // Write WENO results
    std::ofstream weno_file(output_path / config_weno.weno_filename);
    weno_file << "# x, u (WENO)\n";
    weno_file << std::fixed << std::setprecision(6);
    for (int i = 0; i < solver.get_domain().mesh.ncells; ++i) {
        weno_file << std::setw(12) << solver.get_domain().mesh.xcc[i] 
                 << std::setw(12) << u_weno[i] << "\n";
    }
    weno_file.close();
    
    // Write analytical results
    std::ofstream analytical_file(output_path / config_weno.analytical_filename);
    analytical_file << "# x, u (Analytical)\n";
    analytical_file << std::fixed << std::setprecision(6);
    for (int i = 0; i < solver.get_domain().mesh.ncells; ++i) {
        analytical_file << std::setw(12) << solver.get_domain().mesh.xcc[i] 
                       << std::setw(12) << u_analytical[i] << "\n";
    }
    analytical_file.close();
    
    // Print statistics
    auto eno_minmax = std::minmax_element(u_eno.begin(), u_eno.end());
    auto weno_minmax = std::minmax_element(u_weno.begin(), u_weno.end());
    auto anal_minmax = std::minmax_element(u_analytical.begin(), u_analytical.end());
    
    std::cout << "\nSimulation statistics:" << std::endl;
    std::cout << "  ENO min/max: " << *eno_minmax.first << " / " << *eno_minmax.second << std::endl;
    std::cout << "  WENO min/max: " << *weno_minmax.first << " / " << *weno_minmax.second << std::endl;
    std::cout << "  Analytical min/max: " << *anal_minmax.first << " / " << *anal_minmax.second << std::endl;
    
    std::cout << "\nResults written to:" << std::endl;
    std::cout << "  " << (output_path / config_eno.eno_filename).string() << std::endl;
    std::cout << "  " << (output_path / config_weno.weno_filename).string() << std::endl;
    std::cout << "  " << (output_path / config_weno.analytical_filename).string() << std::endl;
    
    std::cout << "\nTo generate the comparison plot, run:" << std::endl;
    std::cout << "  python scripts/postprocess.py" << std::endl;
    std::cout << "==========================================" << std::endl;
}

} // namespace cfd