// ==================== src/flux/FluxFactory.cpp ====================
#include "flux/FluxFactory.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace cfd {

// 创建一个辅助函数用于字符串比较
static std::string to_lower(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}

// 将旧版数字flux_type转换为字符串（保持向后兼容）
static std::string convert_flux_type(int flux_type) {
    switch (flux_type) {
    case 0: return "rusanov";
    case 1: return "engquist-osher";
    case 2: return "lax-wendroff";
    case 3: return "upwind";
    case 4: return "central";
    case 5: return "roe";
    case 6: return "hll";
    default: return "rusanov";
    }
}

std::unique_ptr<FluxCalculator> FluxFactory::create_flux_calculator(int flux_type) {
    std::string flux_name = convert_flux_type(flux_type);
    return create_flux_calculator(flux_name);
}

std::unique_ptr<FluxCalculator> FluxFactory::create_flux_calculator(const std::string& flux_name) {
    std::string name_lower = flux_name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(), ::tolower);
    
    if (name_lower == "rusanov" || name_lower == "llf" || name_lower == "lax-friedrichs") {
        return std::make_unique<RusanovFlux>();
    }
    else if (name_lower == "engquist" || name_lower == "engquist-osher" || name_lower == "eo") {
        return std::make_unique<EngquistOsherFlux>();
    }
    else if (name_lower == "lax-wendroff" || name_lower == "lw") {
        return std::make_unique<LaxWendroffFlux>();
    }
    else if (name_lower == "upwind" || name_lower == "godunov") {
        return std::make_unique<UpwindFlux>();
    }
    else if (name_lower == "central" || name_lower == "centered") {
        return std::make_unique<CentralFlux>();
    }
    else if (name_lower == "roe") {
        return std::make_unique<RoeFlux>();
    }
    else if (name_lower == "hll") {
        return std::make_unique<HLLFlux>();
    }
    else {
        throw std::invalid_argument("Unknown flux name: " + flux_name);
    }
}

std::unique_ptr<FluxCalculator> FluxFactory::create_flux_calculator(int flux_type, Real dt, Real dx) {
    std::string flux_name = convert_flux_type(flux_type);
    return create_flux_calculator(flux_name, dt, dx);
}

std::unique_ptr<FluxCalculator> FluxFactory::create_flux_calculator(const std::string& flux_name, 
    Real dt, Real dx) {
    std::string name_lower = to_lower(flux_name);

    if (name_lower == "lax-wendroff" || name_lower == "lw") {
        return std::make_unique<LaxWendroffFlux>(dt, dx);
    } else {
        return create_flux_calculator(flux_name);
    }
}

std::unique_ptr<FluxCalculator> FluxFactory::create_from_config(const CfdConfig& config) {
    std::string flux_lower = to_lower(config.flux_type);

    if (flux_lower == "lax-wendroff" || flux_lower == "lw") {
        // Lax-Wendroff needs dt and dx
        Real dx = (config.xmax - config.xmin) / config.ncells;
        return std::make_unique<LaxWendroffFlux>(config.dt, dx);
    } else {
        return create_flux_calculator(config.flux_type);
    }
}

std::vector<std::string> FluxFactory::available_fluxes() {
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

std::vector<std::string> FluxFactory::available_fluxes_by_category(const std::string& category) {
    std::vector<std::string> result;
    std::vector<std::string> all_fluxes = {
        "Rusanov", "Engquist-Osher", "Lax-Wendroff", 
        "Upwind", "Central", "Roe", "HLL"
    };
    
    for (const auto& flux_name : all_fluxes) {
        auto flux = create_flux_calculator(flux_name);
        if (flux->category() == category) {
            result.push_back(flux_name);
        }
    }
    
    return result;
}

std::string FluxFactory::get_flux_info(int flux_type) {
    try {
        auto flux = create_flux_calculator(flux_type);
        std::ostringstream oss;
        oss << flux->name() << " (ID: " << flux_type << ")\n";
        oss << "  Category: " << flux->category() << "\n";
        oss << "  Upwind: " << (flux->is_upwind() ? "Yes" : "No") << "\n";
        oss << "  Centered: " << (flux->is_centered() ? "Yes" : "No") << "\n";
        oss << "  Has dissipation: " << (flux->has_dissipation() ? "Yes" : "No") << "\n";
        oss << "  Conservative: " << (flux->is_conservative() ? "Yes" : "No") << "\n";
        oss << "  Entropy stable: " << (flux->is_entropy_stable() ? "Yes" : "No");
        return oss.str();
    } catch (...) {
        return "Unknown flux type: " + std::to_string(flux_type);
    }
}

std::string FluxFactory::get_flux_info(const std::string& flux_name) {
    try {
        auto flux = create_flux_calculator(flux_name);
        return get_flux_info(flux->type_id());
    } catch (...) {
        return "Unknown flux name: " + flux_name;
    }
}

bool FluxFactory::is_available(const std::string& flux_name) {
    std::vector<std::string> available = {
        "rusanov", "engquist", "lax-wendroff", "upwind", 
        "central", "roe", "hll", "llf", "eo", "lw", "godunov"
    };
    
    std::string name_lower = flux_name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(), ::tolower);
    
    return std::find(available.begin(), available.end(), name_lower) != available.end();
}

std::unique_ptr<FluxCalculator> FluxFactory::create_default_flux() {
    return std::make_unique<RusanovFlux>();
}

Vector FluxFactory::compute_flux_comparison(Real uL, Real uR, Real wave_speed) {
    Vector fluxes;
    
    // 手动计算各种通量
    fluxes.push_back(RusanovFlux().compute_flux(uL, uR, wave_speed));
    fluxes.push_back(EngquistOsherFlux().compute_flux(uL, uR, wave_speed));
    fluxes.push_back(LaxWendroffFlux().compute_flux(uL, uR, wave_speed));
    fluxes.push_back(UpwindFlux().compute_flux(uL, uR, wave_speed));
    fluxes.push_back(CentralFlux().compute_flux(uL, uR, wave_speed));
    fluxes.push_back(RoeFlux().compute_flux(uL, uR, wave_speed));
    fluxes.push_back(HLLFlux().compute_flux(uL, uR, wave_speed));
    
    return fluxes;
}

void FluxFactory::print_flux_comparison(Real uL, Real uR, Real wave_speed) {
    std::cout << "Flux Comparison for uL = " << uL << ", uR = " << uR 
              << ", wave_speed = " << wave_speed << ":" << std::endl;
    std::cout << "===================================================" << std::endl;
    
    std::vector<std::unique_ptr<FluxCalculator>> fluxes;
    fluxes.push_back(std::make_unique<RusanovFlux>());
    fluxes.push_back(std::make_unique<EngquistOsherFlux>());
    fluxes.push_back(std::make_unique<LaxWendroffFlux>());
    fluxes.push_back(std::make_unique<UpwindFlux>());
    fluxes.push_back(std::make_unique<CentralFlux>());
    fluxes.push_back(std::make_unique<RoeFlux>());
    fluxes.push_back(std::make_unique<HLLFlux>());
    
    std::cout << std::fixed << std::setprecision(6);
    for (const auto& flux : fluxes) {
        Real flux_value = flux->compute_flux(uL, uR, wave_speed);
        std::cout << "  " << std::setw(25) << std::left << flux->name()
                  << ": " << std::setw(12) << flux_value 
                  << " (type = " << flux->type_id() << ")" << std::endl;
    }
    
    // Exact flux (for linear advection)
    Real exact_flux = wave_speed * 0.5 * (uL + uR);
    std::cout << "  " << std::setw(25) << std::left << "Exact (central)"
              << ": " << std::setw(12) << exact_flux << std::endl;
}

void FluxFactory::compare_all_fluxes(Real uL, Real uR, Real wave_speed) {
    std::cout << "\nComprehensive Flux Comparison" << std::endl;
    std::cout << "==============================" << std::endl;
    
    std::vector<std::unique_ptr<FluxCalculator>> fluxes;
    fluxes.push_back(std::make_unique<RusanovFlux>());
    fluxes.push_back(std::make_unique<EngquistOsherFlux>());
    fluxes.push_back(std::make_unique<LaxWendroffFlux>());
    fluxes.push_back(std::make_unique<UpwindFlux>());
    fluxes.push_back(std::make_unique<CentralFlux>());
    fluxes.push_back(std::make_unique<RoeFlux>());
    fluxes.push_back(std::make_unique<HLLFlux>());
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << std::setw(25) << "Flux Name" 
              << std::setw(15) << "Value" 
              << std::setw(15) << "Category" 
              << std::setw(10) << "Upwind" 
              << std::setw(10) << "Entropy" << std::endl;
    std::cout << std::string(75, '-') << std::endl;
    
    for (const auto& flux : fluxes) {
        Real flux_value = flux->compute_flux(uL, uR, wave_speed);
        std::cout << std::setw(25) << flux->name()
                  << std::setw(15) << flux_value
                  << std::setw(15) << flux->category()
                  << std::setw(10) << (flux->is_upwind() ? "Yes" : "No")
                  << std::setw(10) << (flux->is_entropy_stable() ? "Yes" : "No") << std::endl;
    }
}

} // namespace cfd