#ifndef TIME_INTEGRATOR_FACTORY_H
#define TIME_INTEGRATOR_FACTORY_H

#include "TimeIntegrator.h"
#include "rk1.h"
#include "rk2.h"
#include "rk3.h"
#include <memory>
#include <vector>
#include <string>

namespace cfd {

// ===================================================================
// Time integrator factory
// ===================================================================
class TimeIntegratorFactory {
public:
    // Create integrator by name
    static std::unique_ptr<TimeIntegrator> create_integrator(
        const std::string& method,
        std::unique_ptr<ResidualCalculator> residual_calculator);
    
    // Create integrator by order
    static std::unique_ptr<TimeIntegrator> create_integrator(
        int order,
        std::unique_ptr<ResidualCalculator> residual_calculator);
    
    // Create integrator with configuration
    static std::unique_ptr<TimeIntegrator> create_integrator(
        const CfdConfig& config,
        std::unique_ptr<ResidualCalculator> residual_calculator);
    
    // Get available integrators
    static std::vector<std::string> available_integrators();
    
    // Get integrator info
    static std::string get_integrator_info(const std::string& method);
    
    // Check if integrator is available
    static bool is_available(const std::string& method);
    
    // Get default integrator
    static std::unique_ptr<TimeIntegrator> create_default_integrator(
        std::unique_ptr<ResidualCalculator> residual_calculator);
};

} // namespace cfd

#endif // TIME_INTEGRATOR_FACTORY_H