#ifndef UPWIND_FLUX_H
#define UPWIND_FLUX_H

#include "../FluxCalculator.h"

namespace cfd {

// ===================================================================
// Upwind flux (first-order)
// ===================================================================
class UpwindFlux : public FluxCalculator {
public:
    UpwindFlux() = default;
    ~UpwindFlux() override = default;
    
    Real compute_flux(Real uL, Real uR, Real wave_speed) const override;
    void compute_flux_array(const Vector& uL, const Vector& uR, 
                           Vector& flux, Real wave_speed) const override;
    
    std::string name() const override { return "Upwind"; }
    int type_id() const override { return 3; }
    std::string category() const override { return "Upwind"; }
    bool is_upwind() const override { return true; }
    bool has_dissipation() const override { return false; }
};

} // namespace cfd

#endif // UPWIND_FLUX_H