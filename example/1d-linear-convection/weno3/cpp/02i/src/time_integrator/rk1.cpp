#include "time_integrator/rk1.h"

namespace cfd {

RK1Integrator::RK1Integrator(std::unique_ptr<ResidualCalculator> residual_calculator)
    : TimeIntegrator(std::move(residual_calculator)) {}

void RK1Integrator::step(Vector& u_with_ghosts, 
                       Real dt,
                       const ComputationalDomain& domain) const {
    
    // Get domain parameters
    const int ncells = domain.mesh().ncells();
    const int ist = domain.ist();
    Real dx = domain.mesh().dx();
    
    // Compute residual
    Vector residual(ncells);
    residual_calculator_->compute(u_with_ghosts, residual, domain, dx);
    
    // Update solution (only interior cells)
    for (int i = 0; i < ncells; ++i) {
        u_with_ghosts[ist + i] += dt * residual[i];
    }
}

std::string RK1Integrator::name() const { 
    return "Forward Euler (RK1)";
}

int RK1Integrator::order() const { 
    return 1; 
}

} // namespace cfd