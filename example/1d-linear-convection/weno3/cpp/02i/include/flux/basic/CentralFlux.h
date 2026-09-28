#ifndef CENTRAL_FLUX_H
#define CENTRAL_FLUX_H

#include "../FluxCalculator.h"

namespace cfd {

// ===================================================================
// Central flux
// ===================================================================
class CentralFlux : public FluxCalculator {
public:
    CentralFlux() = default;
    ~CentralFlux() override = default;
    
    Real compute_flux(Real uL, Real uR, Real wave_speed) const override;
    void compute_flux_array(const Vector& uL, const Vector& uR, 
                           Vector& flux, Real wave_speed) const override;
    
    std::string name() const override { return "Central"; }
    int type_id() const override { return 4; }
    std::string category() const override { return "Centered"; }
    bool is_centered() const override { return true; }
    bool has_dissipation() const override { return false; }
    bool is_entropy_stable() const override { return false; } // Central flux is not entropy stable
};

} // namespace cfd

#endif // CENTRAL_FLUX_H