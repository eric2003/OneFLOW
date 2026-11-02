#include "time_integrator/rk2.h"

namespace cfd {

RK2Integrator::RK2Integrator(std::unique_ptr<ResidualCalculator> residual_calculator)
    : TimeIntegrator(std::move(residual_calculator)) {}

void RK2Integrator::step(Vector& u_with_ghosts, 
                       Real dt,
                       const ComputationalDomain& domain) const {
    
    // Get domain parameters
    const int ncells = domain.mesh().ncells();
    const int ist = domain.ist();
    Real dx = domain.mesh().dx();
    
    // Stage 1
    Vector residual(ncells);
    residual_calculator_->compute(u_with_ghosts, residual, domain, dx);
    
    // Temporary storage for stage 1 solution
    Vector u_stage = u_with_ghosts;
    
    // Update to stage 1
    for (int i = 0; i < ncells; ++i) {
        u_stage[ist + i] = u_with_ghosts[ist + i] + dt * residual[i];
    }
    
    // Stage 2
    residual_calculator_->compute(u_stage, residual, domain, dx);
    
    // Final update
    for (int i = 0; i < ncells; ++i) {
        u_with_ghosts[ist + i] = 0.5 * u_with_ghosts[ist + i] + 
                                0.5 * u_stage[ist + i] + 
                                0.5 * dt * residual[i];
    }
}

std::string RK2Integrator::name() const { 
    return "Runge-Kutta 2 (Midpoint)";
}

int RK2Integrator::order() const { 
    return 2; 
}

} // namespace cfd