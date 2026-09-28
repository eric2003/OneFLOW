#ifndef RK3_INTEGRATOR_H
#define RK3_INTEGRATOR_H

#include "TimeIntegrator.h"

namespace cfd {

// ===================================================================
// TVD RK3 time integrator
// ===================================================================
class RK3Integrator : public TimeIntegrator {
public:
    RK3Integrator(std::unique_ptr<ResidualCalculator> residual_calculator);
    ~RK3Integrator() override = default;
    
    void step(Vector& u_with_ghosts, 
             Real dt,
             const ComputationalDomain& domain) const override;
    
    std::string name() const override;
    int order() const override;
    
    // TVD RK3 is TVD preserving
    bool is_tvd() const override { return true; }
    
    // TVD RK3 has CFL limit of 1.0
    Real cfl_limit() const override { return 1.0; }
};

} // namespace cfd

#endif // RK3_INTEGRATOR_H