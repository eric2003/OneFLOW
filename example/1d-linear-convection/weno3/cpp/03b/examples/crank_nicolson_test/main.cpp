// ==================== examples/crank_nicolson_test/main.cpp ====================
#include "cfd_solver.hpp"
#include "analysis.hpp"
#include <iostream>
#include <vector>
// 平台检测
#ifdef _WIN32
#include <windows.h>
#endif

// 跨平台设置控制台编码
void setupConsole() {
#ifdef _WIN32
    // Windows: 设置输入输出为 UTF-8
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);

    // 启用 ANSI 转义序列（Windows 10+）
    HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD mode = 0;
    GetConsoleMode(hOut, &mode);
    SetConsoleMode(hOut, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
#else
    // Linux/macOS: 设置 C 和 C++ 的 locale
    std::setlocale(LC_ALL, "");
    std::locale::global(std::locale(""));
    std::cout.imbue(std::locale());
#endif
}

// Crank-Nicolson测试配置
cfd::CfdConfig create_crank_nicolson_config(int cells = 200, double final_time = 2.0) {
    cfd::CfdConfig config;
    config.ic_type = "complex";
    config.ic_params = {0.5, -0.7, 0.005, 10.0};
    config.recon_scheme = "weno";
    config.spatial_order = 3;
    config.flux_type = "central";  // 中心差分配合C-N
    config.time_integration = "crank-nicolson";
    config.theta = 0.5;
    config.wave_speed = 1.0;
    config.final_time = final_time;
    config.dt = 0.05;  // 更大的时间步长
    config.cfl = 0.5;
    config.ncells = cells;
    config.xmin = -1.0;
    config.xmax = 1.0;
    config.verbose = true;
    return config;
}

cfd::CfdConfig create_rk_config(const std::string& method, 
                               int cells = 200, 
                               double final_time = 2.0) {
    auto config = create_crank_nicolson_config(cells, final_time);
    config.time_integration = method;
    config.flux_type = "rusanov";
    config.dt = 0.001;  // 显式方法需要小时间步长
    return config;
}

int main() {
    setupConsole();

    std::cout << "==========================================" << std::endl;
    std::cout << "Crank-Nicolson Method Test" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // 定义各种方法的配置
    //auto cn_config = create_crank_nicolson_config(200, 2.0);
    auto cn_config = create_crank_nicolson_config(200, 0.5);
    
    // 运行C-N模拟
    std::cout << "\nRunning Crank-Nicolson simulation..." << std::endl;
    auto u_cn = cfd::CfdSolver::run_single_simulation(cn_config);
    
   
    // 计算解析解
    cfd::CfdSolver solver(cn_config);
    auto u_analytical = solver.compute_analytical_solution(cn_config.final_time);
    auto u1 = solver.compute_analytical_solution(0.2);
    
    // 计算误差
    std::cout << "\nError comparison:" << std::endl;
    double cn_error = cfd::Analysis::compute_l1_error(u_cn, u_analytical);
    std::cout << "  Crank-Nicolson: " << cn_error << " (dt=" << cn_config.dt << ")" << std::endl;
    
    // 写入结果
    solver.write_solution("crank_nicolson_results.txt", u_cn);
    solver.write_solution("analytical_results.txt", u_analytical);
    solver.write_solution("u1.txt", u1);
    
    std::cout << "\n✅ Crank-Nicolson test completed!" << std::endl;
    return 0;
}