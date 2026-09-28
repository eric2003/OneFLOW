#include "flux/basic/EngquistOsherFlux.h"
#include <cmath>
#include <stdexcept>

namespace cfd {

Real EngquistOsherFlux::compute_flux(Real uL, Real uR, Real wave_speed) const {
    Real cp = 0.5 * (wave_speed + std::abs(wave_speed));
    Real cm = 0.5 * (wave_speed - std::abs(wave_speed));
    return cp * uL + cm * uR;
}

void EngquistOsherFlux::compute_flux_array(const Vector& uL, const Vector& uR, 
                                         Vector& flux, Real wave_speed) const {
    if (uL.size() != uR.size() || uL.size() != flux.size()) {
        throw std::invalid_argument("Flux array sizes do not match");
    }
    
    Real cp = 0.5 * (wave_speed + std::abs(wave_speed));
    Real cm = 0.5 * (wave_speed - std::abs(wave_speed));
    
    for (size_t i = 0; i < flux.size(); ++i) {
        flux[i] = cp * uL[i] + cm * uR[i];
    }
}

} // namespace cfd