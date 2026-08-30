#include "flux/basic/UpwindFlux.h"
#include <stdexcept>

namespace cfd {

Real UpwindFlux::compute_flux(Real uL, Real uR, Real wave_speed) const {
    if (wave_speed > 0.0) {
        return wave_speed * uL;
    } else {
        return wave_speed * uR;
    }
}

void UpwindFlux::compute_flux_array(const Vector& uL, const Vector& uR, 
                                  Vector& flux, Real wave_speed) const {
    if (uL.size() != uR.size() || uL.size() != flux.size()) {
        throw std::invalid_argument("Flux array sizes do not match");
    }
    
    if (wave_speed > 0.0) {
        for (size_t i = 0; i < flux.size(); ++i) {
            flux[i] = wave_speed * uL[i];
        }
    } else {
        for (size_t i = 0; i < flux.size(); ++i) {
            flux[i] = wave_speed * uR[i];
        }
    }
}

} // namespace cfd