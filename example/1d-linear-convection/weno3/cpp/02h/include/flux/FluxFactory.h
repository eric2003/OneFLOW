#ifndef FLUX_FACTORY_H
#define FLUX_FACTORY_H

#include "config.hpp"
#include "FluxCalculator.h"
#include "basic/RusanovFlux.h"
#include "basic/EngquistOsherFlux.h"
#include "basic/UpwindFlux.h"
#include "basic/CentralFlux.h"
#include "basic/LaxWendroffFlux.h"
#include "advanced/RoeFlux.h"
#include "advanced/HLLFlux.h"
#include <memory>
#include <vector>
#include <string>

namespace cfd {

// ===================================================================
// Flux factory
// ===================================================================
class FluxFactory {
public:
    // Create flux calculator by type ID
    static std::unique_ptr<FluxCalculator> create_flux_calculator(int flux_type);
    
    // Create flux calculator by name
    static std::unique_ptr<FluxCalculator> create_flux_calculator(const std::string& flux_name);
    
    // Create flux calculator with parameters
    static std::unique_ptr<FluxCalculator> create_flux_calculator(int flux_type, 
                                                                 Real dt, Real dx);
    static std::unique_ptr<FluxCalculator> create_flux_calculator(const std::string& flux_name, 
                                                                 Real dt, Real dx);
    
    // Create flux calculator from configuration
    static std::unique_ptr<FluxCalculator> create_from_config(const CfdConfig& config);
    
    // Get available fluxes
    static std::vector<std::string> available_fluxes();
    static std::vector<std::string> available_fluxes_by_category(const std::string& category);
    
    // Get flux information
    static std::string get_flux_info(int flux_type);
    static std::string get_flux_info(const std::string& flux_name);
    
    // Check if flux is available
    static bool is_available(const std::string& flux_name);
    
    // Get default flux
    static std::unique_ptr<FluxCalculator> create_default_flux();
    
    // Flux comparison utilities
    static Vector compute_flux_comparison(Real uL, Real uR, Real wave_speed);
    static void print_flux_comparison(Real uL, Real uR, Real wave_speed);
    static void compare_all_fluxes(Real uL, Real uR, Real wave_speed);
};

} // namespace cfd

#endif // FLUX_FACTORY_H