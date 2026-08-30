#ifndef LAX_WENDROFF_FLUX_H
#define LAX_WENDROFF_FLUX_H

#include "../FluxCalculator.h"

namespace cfd {

// ===================================================================
// Lax-Wendroff flux base class
// ===================================================================
class LaxWendroffFlux : public FluxCalculator {
protected:
    Real dt_ = 0.025;
    Real dx_ = 1.0;
    
public:
    LaxWendroffFlux() = default;
    LaxWendroffFlux(Real dt, Real dx) : dt_(dt), dx_(dx) {}
    ~LaxWendroffFlux() override = default;
    
    Real compute_flux(Real uL, Real uR, Real wave_speed) const override;
    void compute_flux_array(const Vector& uL, const Vector& uR, 
                           Vector& flux, Real wave_speed) const override;
    
    std::string name() const override;
    int type_id() const override { return 2; }
    std::string category() const override { return "Centered"; }
    bool is_centered() const override { return true; }
    bool has_dissipation() const override { return true; }
    
    Real dt() const { return dt_; }
    Real dx() const { return dx_; }
};

} // namespace cfd

#endif // LAX_WENDROFF_FLUX_H