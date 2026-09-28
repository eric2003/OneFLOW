// ==================== src/main.cpp (更新版) ====================
#include "cfd_solver.hpp"
#include "analysis.hpp"
#include <iostream>

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
    cfd::Analysis::run_eno_weno(40, 0.625, true);
    
    std::cout << "\nProgram finished successfully!" << std::endl;
    
    return 0;
}