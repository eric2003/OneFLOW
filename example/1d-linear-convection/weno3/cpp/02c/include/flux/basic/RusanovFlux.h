#ifndef RUSANOV_FLUX_H
#define RUSANOV_FLUX_H

#include "../FluxCalculator.h"

namespace cfd {

// ===================================================================
// Rusanov (Local Lax-Friedrichs) flux
// ===================================================================
class RusanovFlux : public FluxCalculator {
public:
    RusanovFlux() = default;
    ~RusanovFlux() override = default;
    
    Real compute_flux(Real uL, Real uR, Real wave_speed) const override;
    void compute_flux_array(const Vector& uL, const Vector& uR, 
                           Vector& flux, Real wave_speed) const override;
    
    std::string name() const override { return "Rusanov (Local Lax-Friedrichs)"; }
    int type_id() const override { return 0; }
    std::string category() const override { return "Dissipative"; }
    bool has_dissipation() const override { return true; }
    Real dissipation_coefficient() const override { return 1.0; }
};

} // namespace cfd

#endif // RUSANOV_FLUX_H