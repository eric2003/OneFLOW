// ==================== src/analysis.cpp ====================
#include "analysis.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <sstream>

namespace cfd {

void Analysis::run_eno_weno(int cells, double final_time, bool verbose) {
    std::cout << "==========================================" << std::endl;
    std::cout << "OneFLOW-CFD Solver: ENO3 vs WENO3 Analysis" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // Create configurations
    auto [config_eno, config_weno] = CfdConfig::create_analysis_configs(cells, final_time);
    config_eno.output_dir = ".";
    config_weno.output_dir = ".";
    config_eno.verbose = verbose;
    config_weno.verbose = verbose;
    
    // Run ENO simulation
    std::cout << "\nRunning ENO simulation..." << std::endl;
    Vector u_eno = CfdSolver::run_single_simulation(config_eno);
    
    // Run WENO simulation
    std::cout << "\nRunning WENO simulation..." << std::endl;
    Vector u_weno = CfdSolver::run_single_simulation(config_weno);
    
    // Create solver for analytical solution
    CfdSolver solver(config_weno);
    Vector u_analytical = solver.compute_analytical_solution(config_weno.final_time);
    
    // Write results
    std::cout << "\nWriting results to files..." << std::endl;
    
    solver.write_solution("eno_results.txt", u_eno);
    solver.write_solution("weno_results.txt", u_weno);
    solver.write_solution("analytical_results.txt", u_analytical);
    
    // Compute errors
    double eno_error = compute_l1_error(u_eno, u_analytical);
    double weno_error = compute_l1_error(u_weno, u_analytical);
    
    // Print statistics
    auto eno_minmax = std::minmax_element(u_eno.begin(), u_eno.end());
    auto weno_minmax = std::minmax_element(u_weno.begin(), u_weno.end());
    auto anal_minmax = std::minmax_element(u_analytical.begin(), u_analytical.end());
    
    std::cout << "\nSimulation statistics:" << std::endl;
    std::cout << "  ENO min/max: " << *eno_minmax.first << " / " << *eno_minmax.second << std::endl;
    std::cout << "  WENO min/max: " << *weno_minmax.first << " / " << *weno_minmax.second << std::endl;
    std::cout << "  Analytical min/max: " << *anal_minmax.first << " / " << *anal_minmax.second << std::endl;
    
    std::cout << "\nError analysis (L1 norm):" << std::endl;
    std::cout << "  ENO error: " << eno_error << std::endl;
    std::cout << "  WENO error: " << weno_error << std::endl;
    std::cout << "  Error ratio (ENO/WENO): " << eno_error / weno_error << std::endl;
    
    std::cout << "\nResults written to:" << std::endl;
    std::cout << "  eno_results.txt" << std::endl;
    std::cout << "  weno_results.txt" << std::endl;
    std::cout << "  analytical_results.txt" << std::endl;
    
    // Create analysis summary
    std::ofstream summary("analysis_summary.txt");
    summary << "ENO vs WENO Analysis Summary" << std::endl;
    summary << "=============================" << std::endl;
    summary << "Cells: " << cells << std::endl;
    summary << "Final time: " << final_time << std::endl;
    summary << "ENO L1 error: " << eno_error << std::endl;
    summary << "WENO L1 error: " << weno_error << std::endl;
    summary << "Error ratio: " << eno_error / weno_error << std::endl;
    summary.close();
    
    std::cout << "\nAnalysis summary written to: analysis_summary.txt" << std::endl;
    
    std::cout << "\nTo generate the comparison plot, run:" << std::endl;
    std::cout << "  python scripts/postprocess.py" << std::endl;
    std::cout << "==========================================" << std::endl;
}

// 在analysis.cpp的适当位置更新
void Analysis::run_complex_wave_analysis(int cells, double final_time, bool verbose) {
    std::cout << "==========================================" << std::endl;
    std::cout << "Complex Wave Analysis (Precise Parameters)" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // 使用工厂创建复杂波配置
    auto config = cfd::CfdConfig::create_complex_config(cells, final_time);
    config.verbose = verbose;
    
    std::cout << "\nConfiguration:" << std::endl;
    config.print();
    
    // 运行模拟
    std::cout << "\nRunning simulation..." << std::endl;
    cfd::CfdSolver solver(config);
    
    // 现在solver会使用InitialConditionFactory自动创建合适的初始条件
    Vector u_numerical = solver.run_simulation();
    
    // 写入结果
    solver.write_solution("complex_wave_results.txt", u_numerical);
    
    std::cout << "\n✅ Simulation completed!" << std::endl;
    std::cout << "Results written to: complex_wave_results.txt" << std::endl;
    
    // 如果final_time > 0，也可以计算解析解
    if (final_time > 0.0) {
        Vector u_analytical = solver.compute_analytical_solution(final_time);
        solver.write_solution("complex_wave_analytical.txt", u_analytical);
        
        double error = compute_l1_error(u_numerical, u_analytical);
        std::cout << "  L1 error: " << error << std::endl;
    }
}

void Analysis::run_convergence_study(const std::string& scheme,
                                   const std::vector<int>& cell_counts,
                                   double final_time) {
    std::cout << "==========================================" << std::endl;
    std::cout << "Convergence Study: " << scheme << std::endl;
    std::cout << "==========================================" << std::endl;
    
    std::vector<double> errors;
    std::vector<double> resolutions;
    
    for (int cells : cell_counts) {
        std::cout << "\nRunning with " << cells << " cells..." << std::endl;
        
        // Create configuration
        CfdConfig config;
        if (scheme == "eno") {
            config = CfdConfig::create_eno_config(3, cells, final_time);
        } else if (scheme == "weno") {
            config = CfdConfig::create_weno_config(cells, final_time);
        } else {
            throw std::invalid_argument("Unknown scheme: " + scheme);
        }
        
        // Run simulation
        Vector u_numerical = CfdSolver::run_single_simulation(config);
        
        // Compute analytical solution
        CfdSolver solver(config);
        Vector u_analytical = solver.compute_analytical_solution(final_time);
        
        // Compute error
        double error = compute_l1_error(u_numerical, u_analytical);
        errors.push_back(error);
        resolutions.push_back(1.0 / cells);
        
        std::cout << "  L1 error: " << error << std::endl;
    }
    
    // Compute convergence rates
    std::cout << "\nConvergence rates:" << std::endl;
    for (size_t i = 1; i < errors.size(); ++i) {
        double rate = std::log(errors[i-1] / errors[i]) / 
                     std::log(resolutions[i-1] / resolutions[i]);
        std::cout << "  " << cell_counts[i-1] << " → " << cell_counts[i] 
                  << ": " << rate << std::endl;
    }
    
    // Write results
    write_analysis_results("convergence_" + scheme + ".txt", errors, cell_counts);
}

void Analysis::run_flux_comparison(const std::string& recon_scheme,
                                 int cells,
                                 double final_time) {
    std::cout << "==========================================" << std::endl;
    std::cout << "Flux Comparison with " << recon_scheme << std::endl;
    std::cout << "==========================================" << std::endl;
    
    std::vector<std::string> fluxes = {"rusanov", "engquist", "upwind", "central"};
    std::vector<double> errors;
    std::vector<std::string> flux_names;
    
    for (const auto& flux : fluxes) {
        std::cout << "\nTesting " << flux << " flux..." << std::endl;
        
        // Create configuration
        CfdConfig config;
        if (recon_scheme == "eno") {
            config = CfdConfig::create_eno_config(3, cells, final_time);
        } else {
            config = CfdConfig::create_weno_config(cells, final_time);
        }
        
        // Set flux type
        if (flux == "rusanov") config.flux_type = 0;
        else if (flux == "engquist") config.flux_type = 1;
        else if (flux == "upwind") config.flux_type = 3;
        else if (flux == "central") config.flux_type = 4;
        
        // Run simulation
        Vector u_numerical = CfdSolver::run_single_simulation(config);
        
        // Compute analytical solution
        CfdSolver solver(config);
        Vector u_analytical = solver.compute_analytical_solution(final_time);
        
        // Compute error
        double error = compute_l1_error(u_numerical, u_analytical);
        errors.push_back(error);
        flux_names.push_back(flux);
        
        std::cout << "  L1 error: " << error << std::endl;
        
        // Write individual result
        std::string filename = recon_scheme + "_" + flux + "_results.txt";
        solver.write_solution(filename, u_numerical);
    }
    
    // Find best flux
    auto min_iter = std::min_element(errors.begin(), errors.end());
    int best_idx = std::distance(errors.begin(), min_iter);
    
    std::cout << "\nBest performing flux: " << flux_names[best_idx] 
              << " (error = " << errors[best_idx] << ")" << std::endl;
    
    // Write comparison summary
    std::ofstream summary("flux_comparison_" + recon_scheme + ".txt");
    summary << "Flux Comparison for " << recon_scheme << std::endl;
    summary << "==================================" << std::endl;
    for (size_t i = 0; i < fluxes.size(); ++i) {
        summary << std::setw(15) << flux_names[i] << ": " << errors[i] << std::endl;
    }
    summary.close();
}

double Analysis::compute_error(const Vector& numerical, const Vector& analytical) {
    return compute_l1_error(numerical, analytical);
}

double Analysis::compute_l1_error(const Vector& numerical, const Vector& analytical) {
    if (numerical.size() != analytical.size()) {
        throw std::invalid_argument("Vector sizes do not match");
    }
    
    double error = 0.0;
    for (size_t i = 0; i < numerical.size(); ++i) {
        error += std::abs(numerical[i] - analytical[i]);
    }
    return error / numerical.size();
}

double Analysis::compute_l2_error(const Vector& numerical, const Vector& analytical) {
    if (numerical.size() != analytical.size()) {
        throw std::invalid_argument("Vector sizes do not match");
    }
    
    double error = 0.0;
    for (size_t i = 0; i < numerical.size(); ++i) {
        double diff = numerical[i] - analytical[i];
        error += diff * diff;
    }
    return std::sqrt(error / numerical.size());
}

double Analysis::compute_linf_error(const Vector& numerical, const Vector& analytical) {
    if (numerical.size() != analytical.size()) {
        throw std::invalid_argument("Vector sizes do not match");
    }
    
    double max_error = 0.0;
    for (size_t i = 0; i < numerical.size(); ++i) {
        double error = std::abs(numerical[i] - analytical[i]);
        if (error > max_error) max_error = error;
    }
    return max_error;
}

void Analysis::write_analysis_results(const std::string& filename,
                                    const std::vector<double>& errors,
                                    const std::vector<int>& cells) {
    std::ofstream file(filename);
    file << "# Convergence analysis" << std::endl;
    file << "# Cells, Error" << std::endl;
    
    for (size_t i = 0; i < errors.size(); ++i) {
        file << cells[i] << ", " << errors[i] << std::endl;
    }
    
    file.close();
}

void Analysis::plot_results_script(const std::string& output_dir) {
    // Create a Python script for plotting
    std::ofstream script("plot_analysis.py");
    
    script << "import numpy as np\n";
    script << "import matplotlib.pyplot as plt\n";
    script << "import os\n\n";
    
    script << "def read_data(filename):\n";
    script << "    data = np.loadtxt(filename, comments='#')\n";
    script << "    return data[:, 0], data[:, 1]\n\n";
    
    script << "def main():\n";
    script << "    # Read all data files\n";
    script << "    x_eno, u_eno = read_data('eno_results.txt')\n";
    script << "    x_weno, u_weno = read_data('weno_results.txt')\n";
    script << "    x_analytical, u_analytical = read_data('analytical_results.txt')\n\n";
    
    script << "    # Create plot\n";
    script << "    plt.figure(figsize=(12, 8))\n\n";
    
    script << "    # Plot results\n";
    script << "    plt.subplot(2, 1, 1)\n";
    script << "    plt.plot(x_eno, u_eno, 'bo-', linewidth=1, markersize=4, \n";
    script << "             markerfacecolor='none', label='ENO3')\n";
    script << "    plt.plot(x_weno, u_weno, 'gs-', linewidth=1, markersize=4,\n";
    script << "             markerfacecolor='none', label='WENO3')\n";
    script << "    plt.plot(x_analytical, u_analytical, 'r-', linewidth=2, label='Analytical')\n";
    script << "    plt.title('1D Convection: ENO3 vs WENO3')\n";
    script << "    plt.xlabel('x')\n";
    script << "    plt.ylabel('u')\n";
    script << "    plt.legend()\n";
    script << "    plt.grid(True, alpha=0.3)\n\n";
    
    script << "    # Plot errors\n";
    script << "    plt.subplot(2, 1, 2)\n";
    script << "    plt.plot(x_eno, u_eno - u_analytical, 'b-', linewidth=1, label='ENO3 Error')\n";
    script << "    plt.plot(x_weno, u_weno - u_analytical, 'g-', linewidth=1, label='WENO3 Error')\n";
    script << "    plt.title('Error Distribution')\n";
    script << "    plt.xlabel('x')\n";
    script << "    plt.ylabel('Error')\n";
    script << "    plt.legend()\n";
    script << "    plt.grid(True, alpha=0.3)\n\n";
    
    script << "    plt.tight_layout()\n";
    script << "    plt.savefig('" << output_dir << "/analysis_plot.png', dpi=150)\n";
    script << "    plt.show()\n\n";
    
    script << "if __name__ == '__main__':\n";
    script << "    main()\n";
    
    script.close();
    
    std::cout << "Python plotting script created: plot_analysis.py" << std::endl;
}

// 其他方法的实现可以留空或根据需要实现
void Analysis::run_time_integrator_comparison(const std::string& recon_scheme,
                                            int cells,
                                            double final_time) {
    // 实现时间积分器比较
}

void Analysis::run_sensitivity_analysis(const CfdConfig& base_config,
                                      const std::string& parameter,
                                      const std::vector<double>& values) {
    // 实现敏感性分析
}

} // namespace cfd