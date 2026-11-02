#ifndef RK1_INTEGRATOR_H
#define RK1_INTEGRATOR_H

#include "TimeIntegrator.h"

namespace cfd {

// ===================================================================
// Forward Euler (RK1) time integrator
// ===================================================================
class RK1Integrator : public TimeIntegrator {
public:
    RK1Integrator(std::unique_ptr<ResidualCalculator> residual_calculator);
    ~RK1Integrator() override = default;
    
    void step(Vector& u_with_ghosts, 
             Real dt,
             const ComputationalDomain& domain) const override;
    
    std::string name() const override;
    int order() const override;
    
    // Forward Euler has CFL limit of 1.0
    Real cfl_limit() const override { return 1.0; }
};

} // namespace cfd

#endif // RK1_INTEGRATOR_H