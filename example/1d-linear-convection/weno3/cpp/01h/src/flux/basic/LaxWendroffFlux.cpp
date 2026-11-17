#include "flux/basic/LaxWendroffFlux.h"
#include <sstream>

namespace cfd {

Real LaxWendroffFlux::compute_flux(Real uL, Real uR, Real wave_speed) const {
    Real a = wave_speed;
    Real F_L = a * uL;
    Real F_R = a * uR;
    Real u_avg = 0.5 * (uL + uR);
    Real du = uR - uL;
    
    return 0.5 * (F_L + F_R) - 0.5 * (a * dt_ / dx_) * a * du;
}

void LaxWendroffFlux::compute_flux_array(const Vector& uL, const Vector& uR, 
                                       Vector& flux, Real wave_speed) const {
    if (uL.size() != uR.size() || uL.size() != flux.size()) {
        throw std::invalid_argument("Flux array sizes do not match");
    }
    
    Real a = wave_speed;
    Real factor = 0.5 * a * (dt_ / dx_) * a;
    
    for (size_t i = 0; i < flux.size(); ++i) {
        Real F_L = a * uL[i];
        Real F_R = a * uR[i];
        Real du = uR[i] - uL[i];
        flux[i] = 0.5 * (F_L + F_R) - factor * du;
    }
}

std::string LaxWendroffFlux::name() const { 
    std::ostringstream oss;
    oss << "Lax-Wendroff (dt=" << dt_ << ", dx=" << dx_ << ")";
    return oss.str();
}

} // namespace cfd