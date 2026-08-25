// ==================== src/main.cpp ====================
#include "cfd_solver.hpp"

int main() {
    std::cout << "==========================================" << std::endl;
    std::cout << "OneFLOW-CFD Solver for 1D Convection" << std::endl;
    std::cout << "ENO vs WENO Comparison" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    cfd::CfdSolver::perform_eno_weno_analysis();
    
    std::cout << "Program finished successfully!" << std::endl;
    
    return 0;
}