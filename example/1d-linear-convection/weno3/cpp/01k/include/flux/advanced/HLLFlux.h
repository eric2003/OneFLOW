#ifndef HLL_FLUX_H
#define HLL_FLUX_H

#include "../FluxCalculator.h"

namespace cfd {

// ===================================================================
// HLL flux (Harten-Lax-van Leer)
// ===================================================================
class HLLFlux : public FluxCalculator {
public:
    HLLFlux() = default;
    ~HLLFlux() override = default;
    
    Real compute_flux(Real uL, Real uR, Real wave_speed) const override;
    void compute_flux_array(const Vector& uL, const Vector& uR, 
                           Vector& flux, Real wave_speed) const override;
    
    std::string name() const override { return "HLL"; }
    int type_id() const override { return 6; }
    std::string category() const override { return "Approximate Riemann Solver"; }
    bool has_dissipation() const override { return true; }
    bool is_entropy_stable() const override { return true; }
};

} // namespace cfd

#endif // HLL_FLUX_H