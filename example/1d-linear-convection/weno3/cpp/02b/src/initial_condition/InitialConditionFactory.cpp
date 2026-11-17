#include "initial_condition/InitialConditionFactory.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace cfd {

std::unique_ptr<InitialCondition> InitialConditionFactory::create_initial_condition(
    const std::string& ic_type,
    const std::vector<Real>& params) {
    
    std::string type_lower = ic_type;
    std::transform(type_lower.begin(), type_lower.end(), 
                   type_lower.begin(), ::tolower);
    
    if (type_lower == "step" || type_lower == "stepfunction") {
        if (params.size() >= 4) {
            return std::make_unique<StepInitialCondition>(
                params[0], params[1], params[2], params[3]);
        } else if (params.size() >= 2) {
            return std::make_unique<StepInitialCondition>(
                params[0], params[1]);
        } else {
            return std::make_unique<StepInitialCondition>();
        }
    }
    else if (type_lower == "sinusoidal" || type_lower == "sin" || type_lower == "sine") {
        if (params.size() >= 4) {
            return std::make_unique<SinusoidalInitialCondition>(
                params[0], params[1], params[2], params[3]);
        } else if (params.size() >= 3) {
            return std::make_unique<SinusoidalInitialCondition>(
                params[0], params[1], params[2]);
        } else if (params.size() >= 2) {
            return std::make_unique<SinusoidalInitialCondition>(
                params[0], params[1]);
        } else {
            return std::make_unique<SinusoidalInitialCondition>();
        }
    }
    else if (type_lower == "gaussian" || type_lower == "gauss") {
        if (params.size() >= 4) {
            return std::make_unique<GaussianInitialCondition>(
                params[0], params[1], params[2], params[3]);
        } else if (params.size() >= 3) {
            return std::make_unique<GaussianInitialCondition>(
                params[0], params[1], params[2]);
        } else if (params.size() >= 2) {
            return std::make_unique<GaussianInitialCondition>(
                params[0], params[1]);
        } else {
            return std::make_unique<GaussianInitialCondition>();
        }
    }
	else if (type_lower == "complex" || type_lower == "complexwave" || 
			 type_lower == "precise") {
		if (params.size() >= 4) {
			return std::make_unique<ComplexWaveInitialCondition>(
				params[0], params[1], params[2], params[3]
			);
		} else {
			return std::make_unique<ComplexWaveInitialCondition>();
		}
	}	
    else if (type_lower == "custom") {
        // Custom requires a function, which can't be passed through params
        throw std::invalid_argument("Custom initial condition requires a function pointer. "
                                   "Use create_custom() method instead.");
    }
    else {
        // Check custom registry
        auto& registry = get_custom_registry();
        auto it = registry.find(type_lower);
        if (it != registry.end()) {
            return it->second();
        }
        throw std::invalid_argument("Unknown initial condition: " + ic_type);
    }
}

std::unique_ptr<InitialCondition> InitialConditionFactory::create_from_config(
    const CfdConfig& config) {
    // Default to step function for 1D convection
    return create_step(1.0, 2.0, 0.5, 1.0);
}

std::unique_ptr<InitialCondition> InitialConditionFactory::create_step(
    Real low_val, Real high_val, Real start, Real end) {
    return std::make_unique<StepInitialCondition>(low_val, high_val, start, end);
}

std::unique_ptr<InitialCondition> InitialConditionFactory::create_sinusoidal(
    Real amplitude, Real frequency, Real mean, Real phase) {
    return std::make_unique<SinusoidalInitialCondition>(amplitude, frequency, mean, phase);
}

std::unique_ptr<InitialCondition> InitialConditionFactory::create_gaussian(
    Real amplitude, Real center, Real width, Real background) {
    return std::make_unique<GaussianInitialCondition>(amplitude, center, width, background);
}

std::unique_ptr<InitialCondition> InitialConditionFactory::create_custom(
    const std::function<Real(Real)>& func, const std::string& name) {
    return std::make_unique<CustomInitialCondition>(func, name);
}

static std::vector<std::string> available_initial_conditions() {
    return {
        "Step Function",
        "Sinusoidal", 
        "Gaussian",
        "Complex Wave (Precise)",
        "Custom"
    };
}

std::string InitialConditionFactory::get_ic_info(const std::string& ic_type) {
    try {
        auto ic = create_initial_condition(ic_type);
        std::ostringstream oss;
        oss << ic->name() << " (ID: " << ic->type_id() << ")\n";
        oss << "  Type: " << ic_type;
        return oss.str();
    } catch (const std::exception& e) {
        return "Error: " + std::string(e.what());
    }
}

bool InitialConditionFactory::is_available(const std::string& ic_type) {
    std::vector<std::string> available = {
        "step", "sinusoidal", "gaussian", "custom",
        "stepfunction", "sin", "sine", "gauss"
    };
    
    std::string type_lower = ic_type;
    std::transform(type_lower.begin(), type_lower.end(), type_lower.begin(), ::tolower);
    
    if (std::find(available.begin(), available.end(), type_lower) != available.end()) {
        return true;
    }
    
    // Check custom registry
    auto& registry = get_custom_registry();
    return registry.find(type_lower) != registry.end();
}

std::unique_ptr<InitialCondition> InitialConditionFactory::create_default() {
    return create_step();
}

void InitialConditionFactory::register_custom_type(const std::string& name,
                                                  std::function<std::unique_ptr<InitialCondition>()> creator) {
    std::string name_lower = name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(), ::tolower);
    
    get_custom_registry()[name_lower] = creator;
}

std::map<std::string, std::function<std::unique_ptr<InitialCondition>()>>& 
InitialConditionFactory::get_custom_registry() {
    static std::map<std::string, std::function<std::unique_ptr<InitialCondition>()>> registry;
    return registry;
}

} // namespace cfd