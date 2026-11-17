#include "flux/basic/RusanovFlux.h"
#include <cmath>
#include <stdexcept>

namespace cfd {

Real RusanovFlux::compute_flux(Real uL, Real uR, Real wave_speed) const {
    Real F_L = wave_speed * uL;
    Real F_R = wave_speed * uR;
    Real Smax = std::abs(wave_speed);
    return 0.5 * (F_L + F_R) - 0.5 * Smax * (uR - uL);
}

void RusanovFlux::compute_flux_array(const Vector& uL, const Vector& uR, 
                                   Vector& flux, Real wave_speed) const {
    if (uL.size() != uR.size() || uL.size() != flux.size()) {
        throw std::invalid_argument("Flux array sizes do not match");
    }
    
    Real Smax = std::abs(wave_speed);
    
    for (size_t i = 0; i < flux.size(); ++i) {
        Real F_L = wave_speed * uL[i];
        Real F_R = wave_speed * uR[i];
        flux[i] = 0.5 * (F_L + F_R) - 0.5 * Smax * (uR[i] - uL[i]);
    }
}

} // namespace cfd