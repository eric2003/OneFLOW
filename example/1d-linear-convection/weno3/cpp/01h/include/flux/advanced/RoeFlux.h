#ifndef ROE_FLUX_H
#define ROE_FLUX_H

#include "../FluxCalculator.h"

namespace cfd {

// ===================================================================
// Roe flux (approximate Riemann solver)
// ===================================================================
class RoeFlux : public FluxCalculator {
public:
    RoeFlux() = default;
    ~RoeFlux() override = default;
    
    Real compute_flux(Real uL, Real uR, Real wave_speed) const override;
    void compute_flux_array(const Vector& uL, const Vector& uR, 
                           Vector& flux, Real wave_speed) const override;
    
    std::string name() const override { return "Roe"; }
    int type_id() const override { return 5; }
    std::string category() const override { return "Approximate Riemann Solver"; }
    bool is_upwind() const override { return true; }
    bool is_entropy_stable() const override { return false; } // Roe flux can violate entropy condition
};

} // namespace cfd

#endif // ROE_FLUX_H