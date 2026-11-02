#include "flux/advanced/RoeFlux.h"
#include <stdexcept>

namespace cfd {

Real RoeFlux::compute_flux(Real uL, Real uR, Real wave_speed) const {
    Real a = wave_speed;
    Real F_L = a * uL;
    Real F_R = a * uR;
    
    // For linear advection, Roe average speed is just the wave speed
    Real a_roe = a;
    
    if (a_roe >= 0.0) {
        return F_L;
    } else {
        return F_R;
    }
}

void RoeFlux::compute_flux_array(const Vector& uL, const Vector& uR, 
                               Vector& flux, Real wave_speed) const {
    if (uL.size() != uR.size() || uL.size() != flux.size()) {
        throw std::invalid_argument("Flux array sizes do not match");
    }
    
    Real a = wave_speed;
    
    if (a >= 0.0) {
        for (size_t i = 0; i < flux.size(); ++i) {
            flux[i] = a * uL[i];
        }
    } else {
        for (size_t i = 0; i < flux.size(); ++i) {
            flux[i] = a * uR[i];
        }
    }
}

} // namespace cfd