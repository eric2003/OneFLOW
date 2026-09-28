#ifndef INITIAL_CONDITION_FACTORY_H
#define INITIAL_CONDITION_FACTORY_H

#include "InitialCondition.h"
#include "StepInitialCondition.h"
#include "SinusoidalInitialCondition.h"
#include "GaussianInitialCondition.h"
#include "CustomInitialCondition.h"
#include "config.hpp"
#include <memory>
#include <vector>
#include <string>
#include <map>

namespace cfd {

// ===================================================================
// Initial condition factory
// ===================================================================
class InitialConditionFactory {
public:
    // Create initial condition by type name
    static std::unique_ptr<InitialCondition> create_initial_condition(
        const std::string& ic_type,
        const std::vector<Real>& params = {});
    
    // Create initial condition from configuration
    static std::unique_ptr<InitialCondition> create_from_config(
        const CfdConfig& config);
    
    // Create step function (most common for CFD tests)
    static std::unique_ptr<InitialCondition> create_step(
        Real low_val = 1.0, Real high_val = 2.0,
        Real start = 0.5, Real end = 1.0);
    
    // Create sinusoidal
    static std::unique_ptr<InitialCondition> create_sinusoidal(
        Real amplitude = 0.5, Real frequency = 1.0,
        Real mean = 1.0, Real phase = 0.0);
    
    // Create Gaussian
    static std::unique_ptr<InitialCondition> create_gaussian(
        Real amplitude = 1.0, Real center = 1.0,
        Real width = 0.1, Real background = 1.0);
    
    // Create custom
    static std::unique_ptr<InitialCondition> create_custom(
        const std::function<Real(Real)>& func,
        const std::string& name = "Custom");
    
    // Get available initial conditions
    static std::vector<std::string> available_initial_conditions();
    
    // Get information about an initial condition type
    static std::string get_ic_info(const std::string& ic_type);
    
    // Check if initial condition type is available
    static bool is_available(const std::string& ic_type);
    
    // Create default initial condition (step function)
    static std::unique_ptr<InitialCondition> create_default();
    
    // Register custom initial condition type
    static void register_custom_type(const std::string& name,
                                   std::function<std::unique_ptr<InitialCondition>()> creator);
    
private:
    // Registry for custom initial conditions
    static std::map<std::string, std::function<std::unique_ptr<InitialCondition>()>>& 
    get_custom_registry();
};

} // namespace cfd

#endif // INITIAL_CONDITION_FACTORY_H