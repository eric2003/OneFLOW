// ==================== examples/eno_weno_comparison/main.cpp ====================
#include "cfd_solver.hpp"
#include "analysis.hpp"
#include <iostream>

namespace cfd {

using Real = double;
using Vector = std::vector<Real>;

// 具体算例的配置定义
CfdConfig create_eno_config( int cells = 200, double final_time = 2.0 ) {
    cfd::CfdConfig config;
    config.ic_type = "complex";
    config.ic_params = { 0.5, -0.7, 0.005, 10.0 };
    config.recon_scheme = "eno";
    config.spatial_order = 3;
    config.flux_type = "rusanov";
    config.time_integration = "rk3";
    config.rk_order = 3;
    config.wave_speed = 1.0;
    config.final_time = final_time;
    config.dt = 0.025;
    config.cfl = 0.5;
    config.ncells = cells;
    config.xmin = -1.0;
    config.xmax = 1.0;
    config.verbose = true;
    return config;
}

CfdConfig create_weno_config( int cells = 200, double final_time = 2.0 ) {
    auto config = create_eno_config( cells, final_time );
    config.recon_scheme = "weno";
    return config;
}

void run_eno_weno( int cells, int spatial_order, double final_time, bool verbose ) {
    std::cout << "==========================================" << std::endl;

    // 根据 spatial_order 动态生成标题
    std::string eno_name, weno_name;
    if ( spatial_order == 3 ) {
        eno_name = "ENO3";
        weno_name = "WENO3";
    }
    else if ( spatial_order == 5 ) {
        eno_name = "ENO5";
        weno_name = "WENO5";
    }
    else {
        // 默认使用输入的数字
        eno_name = "ENO" + std::to_string( spatial_order );
        weno_name = "WENO" + std::to_string( spatial_order );
    }

    std::cout << "OneFLOW-CFD Solver: " << eno_name << " vs " << weno_name << " Analysis" << std::endl;
    std::cout << "==========================================" << std::endl;

    // 判断使用哪种初始条件
    std::string ic_type = "complex";  // 这里改为 "complex" 使用复杂波

    // 根据初始条件类型设置合适的网格范围
    double xmin, xmax;
    std::vector<Real> ic_params;

    if ( ic_type == "complex" || ic_type == "complexwave" || ic_type == "precise" ) {
        // 复杂波需要 [-1.0, 1.0] 域
        xmin = -1.0;
        xmax = 1.0;
        // 复杂波精确参数：a, z, delta, alpha
        ic_params = { 0.5, -0.7, 0.005, 10.0 };
        std::cout << "Using Complex Wave initial condition" << std::endl;
    }
    else {
        // 默认使用 step 函数，域为 [0.0, 2.0]
        xmin = 0.0;
        xmax = 2.0;
        // step 函数参数：low_val, high_val, start, end
        ic_params = { 1.0, 2.0, 0.5, 1.0 };
        std::cout << "Using Step Function initial condition" << std::endl;
    }

    // 创建配置
    CfdConfig config_eno;
    config_eno.ic_type = ic_type;
    config_eno.ic_params = ic_params;
    config_eno.recon_scheme = "eno";
    config_eno.spatial_order = spatial_order;  // 使用输入的阶数
    config_eno.flux_type = "rusanov";
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
    Vector u_eno = CfdSolver::run_single_simulation( config_eno );

    // Run WENO simulation
    std::cout << "\nRunning " << weno_name << " simulation..." << std::endl;
    Vector u_weno = CfdSolver::run_single_simulation( config_weno );

    // Create solver for analytical solution
    CfdSolver solver( config_weno );
    Vector u_analytical = solver.compute_analytical_solution( config_weno.final_time );

    // Write results - 根据初始条件类型和阶数使用不同的文件名
    std::string suffix = "_order" + std::to_string( spatial_order );

    std::cout << "\nWriting results to files..." << std::endl;

    solver.write_solution( "eno_results" + suffix + ".txt", u_eno );
    solver.write_solution( "weno_results" + suffix + ".txt", u_weno );
    solver.write_solution( "analytical_results" + suffix + ".txt", u_analytical );

    // Compute errors
    double eno_error = Analysis::compute_l1_error( u_eno, u_analytical );
    double weno_error = Analysis::compute_l1_error( u_weno, u_analytical );

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
    std::ofstream summary( summary_filename );
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
}

int main() {
    std::cout << "==========================================" << std::endl;
    std::cout << "OneFLOW-CFD Solver for 1D Convection" << std::endl;
    std::cout << "==========================================" << std::endl;

    std::cout << "\nAvailable options:" << std::endl;
    std::cout << "  1. Run ENO vs WENO comparison (default)" << std::endl;
    std::cout << "  2. Run convergence study" << std::endl;
    std::cout << "  3. Run flux comparison" << std::endl;
    std::cout << "  4. Run custom simulation" << std::endl;

    // 默认运行ENO vs WENO比较
    //cfd::run_eno_weno(40, 0.625, true);
    cfd::run_eno_weno(200, 5,  8.0, true);

    std::cout << "\nProgram finished successfully!" << std::endl;
    return 0;
}