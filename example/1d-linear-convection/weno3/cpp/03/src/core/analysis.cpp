// ==================== src/analysis.cpp ====================
#include "analysis.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <sstream>

namespace cfd {

void Analysis::run_eno_weno(int cells, int spatial_order, double final_time, bool verbose) {
    std::cout << "==========================================" << std::endl;

    // 根据 spatial_order 动态生成标题
    std::string eno_name, weno_name;
    if (spatial_order == 3) {
        eno_name = "ENO3";
        weno_name = "WENO3";
    } else if (spatial_order == 5) {
        eno_name = "ENO5";
        weno_name = "WENO5";
    } else {
        // 默认使用输入的数字
        eno_name = "ENO" + std::to_string(spatial_order);
        weno_name = "WENO" + std::to_string(spatial_order);
    }

    std::cout << "OneFLOW-CFD Solver: " << eno_name << " vs " << weno_name << " Analysis" << std::endl;
    std::cout << "==========================================" << std::endl;

    // 判断使用哪种初始条件
    std::string ic_type = "complex";  // 这里改为 "complex" 使用复杂波

    // 根据初始条件类型设置合适的网格范围
    double xmin, xmax;
    std::vector<Real> ic_params;

    if (ic_type == "complex" || ic_type == "complexwave" || ic_type == "precise") {
        // 复杂波需要 [-1.0, 1.0] 域
        xmin = -1.0;
        xmax = 1.0;
        // 复杂波精确参数：a, z, delta, alpha
        ic_params = {0.5, -0.7, 0.005, 10.0};
        std::cout << "Using Complex Wave initial condition" << std::endl;
    } else {
        // 默认使用 step 函数，域为 [0.0, 2.0]
        xmin = 0.0;
        xmax = 2.0;
        // step 函数参数：low_val, high_val, start, end
        ic_params = {1.0, 2.0, 0.5, 1.0};
        std::cout << "Using Step Function initial condition" << std::endl;
    }

    // 创建配置
    CfdConfig config_eno;
    config_eno.ic_type = ic_type;
    config_eno.ic_params = ic_params;
    config_eno.recon_scheme = "eno";
    config_eno.spatial_order = spatial_order;  // 使用输入的阶数
    config_eno.flux_type = 0;
    config_eno.rk_order = 3;
    config_eno.wave_speed = 1.0;
    config_eno.final_time = final_time;
    config_eno.dt = 0.025;
    config_eno.cfl = 0.5;
    config_eno.ncells = cells;
    config_eno.xmin = xmin;    // 设置网格范围
    config_eno.xmax = xmax;
    config_eno.output_dir = ".";
    config_eno.verbose = verbose;

    CfdConfig config_weno = config_eno;  // 复制配置
    config_weno.recon_scheme = "weno";   // 修改重构方案

    // 打印配置信息
    std::cout << "\nDomain: [" << xmin << ", " << xmax << "]" << std::endl;
    std::cout << "Cells: " << cells << std::endl;
    std::cout << "Spatial order: " << spatial_order << std::endl;

    // Run ENO simulation
    std::cout << "\nRunning " << eno_name << " simulation..." << std::endl;
    Vector u_eno = CfdSolver::run_single_simulation(config_eno);

    // Run WENO simulation
    std::cout << "\nRunning " << weno_name << " simulation..." << std::endl;
    Vector u_weno = CfdSolver::run_single_simulation(config_weno);

    // Create solver for analytical solution
    CfdSolver solver(config_weno);
    Vector u_analytical = solver.compute_analytical_solution(config_weno.final_time);

    // Write results - 根据初始条件类型和阶数使用不同的文件名
    std::string suffix = "_order" + std::to_string(spatial_order);

    std::cout << "\nWriting results to files..." << std::endl;

    solver.write_solution("eno_results" + suffix + ".txt", u_eno);
    solver.write_solution("weno_results" + suffix + ".txt", u_weno);
    solver.write_solution("analytical_results" + suffix + ".txt", u_analytical);

    // Compute errors
    double eno_error = compute_l1_error(u_eno, u_analytical);
    double weno_error = compute_l1_error(u_weno, u_analytical);

    // Print statistics
    std::cout << "\nError analysis (L1 norm):" << std::endl;
    std::cout << "  " << eno_name << " error: " << eno_error << std::endl;
    std::cout << "  " << weno_name << " error: " << weno_error << std::endl;
    std::cout << "  Error ratio (" << eno_name << "/" << weno_name << "): " << eno_error / weno_error << std::endl;

    std::cout << "\nResults written to:" << std::endl;
    std::cout << "  eno_results" + suffix + ".txt" << std::endl;
    std::cout << "  weno_results" + suffix + ".txt" << std::endl;
    std::cout << "  analytical_results" + suffix + ".txt" << std::endl;

    // 创建分析摘要
    std::string summary_filename = "analysis_summary" + suffix + ".txt";
    std::ofstream summary(summary_filename);
    summary << ic_type << " " << eno_name << " vs " << weno_name << " Analysis Summary" << std::endl;
    summary << "==========================================" << std::endl;
    summary << "Initial condition: " << ic_type << std::endl;
    summary << "Cells: " << cells << std::endl;
    summary << "Domain: [" << xmin << ", " << xmax << "]" << std::endl;
    summary << "Final time: " << final_time << std::endl;
    summary << "Spatial order: " << spatial_order << std::endl;
    summary << eno_name << " L1 error: " << eno_error << std::endl;
    summary << weno_name << " L1 error: " << weno_error << std::endl;
    summary << "Error ratio (" << eno_name << "/" << weno_name << "): " << eno_error / weno_error << std::endl;
    summary.close();

    std::cout << "\nAnalysis summary written to: " << summary_filename << std::endl;

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