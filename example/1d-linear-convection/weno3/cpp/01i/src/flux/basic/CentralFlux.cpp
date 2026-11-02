#include "flux/basic/CentralFlux.h"
#include <stdexcept>

namespace cfd {

Real CentralFlux::compute_flux(Real uL, Real uR, Real wave_speed) const {
    return 0.5 * wave_speed * (uL + uR);
}

void CentralFlux::compute_flux_array(const Vector& uL, const Vector& uR, 
                                   Vector& flux, Real wave_speed) const {
    if (uL.size() != uR.size() || uL.size() != flux.size()) {
        throw std::invalid_argument("Flux array sizes do not match");
    }
    
    for (size_t i = 0; i < flux.size(); ++i) {
        flux[i] = 0.5 * wave_speed * (uL[i] + uR[i]);
    }
}

} // namespace cfd