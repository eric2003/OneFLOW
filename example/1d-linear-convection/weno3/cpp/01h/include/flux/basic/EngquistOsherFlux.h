#ifndef ENGQUIST_OSHER_FLUX_H
#define ENGQUIST_OSHER_FLUX_H

#include "../FluxCalculator.h"

namespace cfd {

// ===================================================================
// Engquist-Osher flux
// ===================================================================
class EngquistOsherFlux : public FluxCalculator {
public:
    EngquistOsherFlux() = default;
    ~EngquistOsherFlux() override = default;
    
    Real compute_flux(Real uL, Real uR, Real wave_speed) const override;
    void compute_flux_array(const Vector& uL, const Vector& uR, 
                           Vector& flux, Real wave_speed) const override;
    
    std::string name() const override { return "Engquist-Osher"; }
    int type_id() const override { return 1; }
    std::string category() const override { return "Upwind"; }
    bool is_upwind() const override { return true; }
    bool has_dissipation() const override { return false; }
    bool is_entropy_stable() const override { return true; }
};

} // namespace cfd

#endif // ENGQUIST_OSHER_FLUX_H