#ifndef TIME_INTEGRATOR_H
#define TIME_INTEGRATOR_H

#include "residual.hpp"
#include "domain.hpp"
#include <memory>
#include <vector>
#include <string>

namespace cfd {

// Type aliases
using Real = double;
using Vector = std::vector<Real>;

// ===================================================================
// Time integrator base class
// ===================================================================
class TimeIntegrator {
protected:
    std::unique_ptr<ResidualCalculator> residual_calculator_;
    
public:
    TimeIntegrator(std::unique_ptr<ResidualCalculator> residual_calculator);
    virtual ~TimeIntegrator() = default;
    
    // Perform one time step
    virtual void step(Vector& u_with_ghosts, 
                     Real dt,
                     const ComputationalDomain& domain) const = 0;
    
    // Get time integrator name
    virtual std::string name() const = 0;
    
    // Get order of accuracy
    virtual int order() const = 0;
    
    // Get residual calculator
    const ResidualCalculator& residual_calculator() const;
    
    // Check if integrator is explicit
    virtual bool is_explicit() const { return true; }
    
    // Check if integrator is implicit
    virtual bool is_implicit() const { return false; }
    
    // Check if integrator is TVD (Total Variation Diminishing)
    virtual bool is_tvd() const { return false; }
    
    // Get stability region (CFL limit for explicit methods)
    virtual Real cfl_limit() const { return 1.0; }
};

} // namespace cfd

#endif // TIME_INTEGRATOR_H