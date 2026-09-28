#include "flux/advanced/HLLFlux.h"
#include <stdexcept>

namespace cfd {

Real HLLFlux::compute_flux(Real uL, Real uR, Real wave_speed) const {
    Real a = wave_speed;
    Real F_L = a * uL;
    Real F_R = a * uR;
    
    // Wave speeds for linear advection
    Real S_L = a;
    Real S_R = a;
    
    if (S_L >= 0.0) {
        return F_L;
    } else if (S_R <= 0.0) {
        return F_R;
    } else {
        return (S_R * F_L - S_L * F_R + S_L * S_R * (uR - uL)) / (S_R - S_L);
    }
}

void HLLFlux::compute_flux_array(const Vector& uL, const Vector& uR, 
                               Vector& flux, Real wave_speed) const {
    if (uL.size() != uR.size() || uL.size() != flux.size()) {
        throw std::invalid_argument("Flux array sizes do not match");
    }
    
    Real a = wave_speed;
    Real F_L, F_R, S_L, S_R;
    
    for (size_t i = 0; i < flux.size(); ++i) {
        F_L = a * uL[i];
        F_R = a * uR[i];
        S_L = a;
        S_R = a;
        
        if (S_L >= 0.0) {
            flux[i] = F_L;
        } else if (S_R <= 0.0) {
            flux[i] = F_R;
        } else {
            flux[i] = (S_R * F_L - S_L * F_R + S_L * S_R * (uR[i] - uL[i])) / (S_R - S_L);
        }
    }
}

} // namespace cfd