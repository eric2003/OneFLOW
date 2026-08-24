// ==================== src/flux.cpp (修复版) ====================
#include "flux.hpp"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>

namespace cfd {

// ===================================================================
// Basic flux implementations
// ===================================================================

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

// LaxWendroffFlux implementation (base class)
Real LaxWendroffFlux::compute_flux(Real uL, Real uR, Real wave_speed) const {
    throw std::runtime_error("LaxWendroffFlux needs dt and dx parameters");
}

void LaxWendroffFlux::compute_flux_array(const Vector& uL, const Vector& uR, 
                                       Vector& flux, Real wave_speed) const {
    throw std::runtime_error("LaxWendroffFlux needs dt and dx parameters");
}

// LaxWendroffFluxExtended implementation
Real LaxWendroffFluxExtended::compute_flux(Real uL, Real uR, Real wave_speed) const {
    Real a = wave_speed;
    Real F_L = a * uL;
    Real F_R = a * uR;
    Real u_avg = 0.5 * (uL + uR);
    Real du = uR - uL;
    
    return 0.5 * (F_L + F_R) - 0.5 * (a * dt_ / dx_) * a * du;
}

void LaxWendroffFluxExtended::compute_flux_array(const Vector& uL, const Vector& uR, 
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

std::string LaxWendroffFluxExtended::name() const { 
    return "Lax-Wendroff (dt=" + std::to_string(dt_) + ", dx=" + std::to_string(dx_) + ")";
}

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

// ===================================================================
// FluxFactory implementation
// ===================================================================

// Basic factory methods
std::unique_ptr<FluxCalculator> FluxFactory::create_flux_calculator(int flux_type) {
    switch (flux_type) {
        case 0:
            return std::make_unique<RusanovFlux>();
        case 1:
            return std::make_unique<EngquistOsherFlux>();
        case 2:
            return std::make_unique<LaxWendroffFlux>();
        case 3:
            return std::make_unique<UpwindFlux>();
        case 4:
            return std::make_unique<CentralFlux>();
        default:
            throw std::invalid_argument("Unknown flux type: " + std::to_string(flux_type));
    }
}

std::unique_ptr<FluxCalculator> FluxFactory::create_flux_calculator(const std::string& flux_name) {
    if (flux_name == "rusanov" || flux_name == "Rusanov") {
        return std::make_unique<RusanovFlux>();
    } else if (flux_name == "engquist" || flux_name == "Engquist-Osher") {
        return std::make_unique<EngquistOsherFlux>();
    } else if (flux_name == "lax-wendroff" || flux_name == "Lax-Wendroff") {
        return std::make_unique<LaxWendroffFlux>();
    } else if (flux_name == "upwind" || flux_name == "Upwind") {
        return std::make_unique<UpwindFlux>();
    } else if (flux_name == "central" || flux_name == "Central") {
        return std::make_unique<CentralFlux>();
    } else {
        throw std::invalid_argument("Unknown flux name: " + flux_name);
    }
}

// Extended factory methods with dt and dx parameters
std::unique_ptr<FluxCalculator> FluxFactory::create_flux_calculator(int flux_type, Real dt, Real dx) {
    switch (flux_type) {
        case 0:
            return std::make_unique<RusanovFlux>();
        case 1:
            return std::make_unique<EngquistOsherFlux>();
        case 2:
            return std::make_unique<LaxWendroffFluxExtended>(dt, dx);
        case 3:
            return std::make_unique<UpwindFlux>();
        case 4:
            return std::make_unique<CentralFlux>();
        default:
            throw std::invalid_argument("Unknown flux type: " + std::to_string(flux_type));
    }
}

std::unique_ptr<FluxCalculator> FluxFactory::create_flux_calculator(const std::string& flux_name, Real dt, Real dx) {
    if (flux_name == "rusanov" || flux_name == "Rusanov") {
        return std::make_unique<RusanovFlux>();
    } else if (flux_name == "engquist" || flux_name == "Engquist-Osher") {
        return std::make_unique<EngquistOsherFlux>();
    } else if (flux_name == "lax-wendroff" || flux_name == "Lax-Wendroff") {
        return std::make_unique<LaxWendroffFluxExtended>(dt, dx);
    } else if (flux_name == "upwind" || flux_name == "Upwind") {
        return std::make_unique<UpwindFlux>();
    } else if (flux_name == "central" || flux_name == "Central") {
        return std::make_unique<CentralFlux>();
    } else {
        throw std::invalid_argument("Unknown flux name: " + flux_name);
    }
}

// Extended factory with all flux types
std::unique_ptr<FluxCalculator> FluxFactory::create_extended_flux_calculator(int flux_type, Real dt, Real dx) {
    switch (flux_type) {
        case 0:
            return std::make_unique<RusanovFlux>();
        case 1:
            return std::make_unique<EngquistOsherFlux>();
        case 2:
            return std::make_unique<LaxWendroffFluxExtended>(dt, dx);
        case 3:
            return std::make_unique<UpwindFlux>();
        case 4:
            return std::make_unique<CentralFlux>();
        case 5:
            return std::make_unique<RoeFlux>();
        case 6:
            return std::make_unique<HLLFlux>();
        default:
            throw std::invalid_argument("Unknown flux type: " + std::to_string(flux_type));
    }
}

// Get available fluxes
std::vector<std::string> FluxFactory::available_fluxes() {
    return {
        "Rusanov (Local Lax-Friedrichs)",
        "Engquist-Osher",
        "Lax-Wendroff",
        "Upwind",
        "Central"
    };
}

std::vector<std::string> FluxFactory::available_fluxes_extended() {
    return {
        "Rusanov (Local Lax-Friedrichs)",
        "Engquist-Osher",
        "Lax-Wendroff",
        "Upwind",
        "Central",
        "Roe",
        "HLL"
    };
}

// ===================================================================
// Flux comparison utility functions (修复版)
// ===================================================================

Vector FluxFactory::compute_flux_comparison(Real uL, Real uR, Real wave_speed) {
    Vector fluxes;
    
    // 手动创建和计算通量，避免 unique_ptr 复制问题
    fluxes.push_back(RusanovFlux().compute_flux(uL, uR, wave_speed));
    fluxes.push_back(EngquistOsherFlux().compute_flux(uL, uR, wave_speed));
    fluxes.push_back(UpwindFlux().compute_flux(uL, uR, wave_speed));
    fluxes.push_back(CentralFlux().compute_flux(uL, uR, wave_speed));
    fluxes.push_back(RoeFlux().compute_flux(uL, uR, wave_speed));
    fluxes.push_back(HLLFlux().compute_flux(uL, uR, wave_speed));
    
    return fluxes;
}

void FluxFactory::print_flux_comparison(Real uL, Real uR, Real wave_speed) {
    std::cout << "Flux Comparison for uL = " << uL << ", uR = " << uR 
              << ", wave_speed = " << wave_speed << ":" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    
    // 使用指针数组，避免 vector 中的 unique_ptr
    struct FluxInfo {
        int type_id;
        std::string name;
        Real flux_value;
    };
    
    std::vector<FluxInfo> flux_infos;
    
    // 计算各种通量
    flux_infos.push_back({0, "Rusanov", RusanovFlux().compute_flux(uL, uR, wave_speed)});
    flux_infos.push_back({1, "Engquist-Osher", EngquistOsherFlux().compute_flux(uL, uR, wave_speed)});
    flux_infos.push_back({3, "Upwind", UpwindFlux().compute_flux(uL, uR, wave_speed)});
    flux_infos.push_back({4, "Central", CentralFlux().compute_flux(uL, uR, wave_speed)});
    flux_infos.push_back({5, "Roe", RoeFlux().compute_flux(uL, uR, wave_speed)});
    flux_infos.push_back({6, "HLL", HLLFlux().compute_flux(uL, uR, wave_speed)});
    
    for (const auto& info : flux_infos) {
        std::cout << "  " << std::setw(25) << std::left << info.name
                  << ": " << std::setw(10) << std::right << info.flux_value 
                  << " (type = " << info.type_id << ")" << std::endl;
    }
    
    // Exact flux (for linear advection)
    Real exact_flux = wave_speed * 0.5 * (uL + uR);  // Central difference
    std::cout << "  " << std::setw(25) << std::left << "Exact (central)"
              << ": " << std::setw(10) << std::right << exact_flux << std::endl;
}

} // namespace cfd