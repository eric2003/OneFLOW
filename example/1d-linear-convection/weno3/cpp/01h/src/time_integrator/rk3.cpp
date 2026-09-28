#include "time_integrator/rk3.h"

namespace cfd {

RK3Integrator::RK3Integrator(std::unique_ptr<ResidualCalculator> residual_calculator)
    : TimeIntegrator(std::move(residual_calculator)) {}

void RK3Integrator::step(Vector& u_with_ghosts, 
                       Real dt,
                       const ComputationalDomain& domain) const {
    
    // Get domain parameters
    const int ncells = domain.mesh().ncells();
    const int ist = domain.ist();
    Real dx = domain.mesh().dx();
    
    Vector residual(ncells);
    
    // Stage 1
    residual_calculator_->compute(u_with_ghosts, residual, domain, dx);
    Vector u1 = u_with_ghosts;
    for (int i = 0; i < ncells; ++i) {
        u1[ist + i] = u_with_ghosts[ist + i] + dt * residual[i];
    }
    
    // Stage 2
    residual_calculator_->compute(u1, residual, domain, dx);
    Vector u2 = u_with_ghosts;
    for (int i = 0; i < ncells; ++i) {
        u2[ist + i] = 0.75 * u_with_ghosts[ist + i] + 
                     0.25 * u1[ist + i] + 
                     0.25 * dt * residual[i];
    }
    
    // Stage 3
    residual_calculator_->compute(u2, residual, domain, dx);
    for (int i = 0; i < ncells; ++i) {
        u_with_ghosts[ist + i] = (1.0/3.0) * u_with_ghosts[ist + i] + 
                                (2.0/3.0) * u2[ist + i] + 
                                (2.0/3.0) * dt * residual[i];
    }
}

std::string RK3Integrator::name() const { 
    return "TVD Runge-Kutta 3";
}

int RK3Integrator::order() const { 
    return 3; 
}

} // namespace cfd