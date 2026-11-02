#include "initial_condition/InitialConditionFactory.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace cfd {

std::unique_ptr<InitialCondition> InitialConditionFactory::create_initial_condition(
    const std::string& ic_type,
    const std::vector<Real>& params) {
    
    std::string type_lower = to_lower(ic_type);
    
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
    else if (type_lower == "complex" || type_lower == "complexwave" || type_lower == "precise") {
        // Complex wave with precise parameters
        if (params.size() >= 4) {
            return std::make_unique<ComplexWaveInitialCondition>(
                params[0], params[1], params[2], params[3]
            );
        } else {
            // Use default precise parameters
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

// ===================================================================
// 从配置创建初始条件（新增）
// ===================================================================

std::unique_ptr<InitialCondition> InitialConditionFactory::create_from_config(
    const CfdConfig& config) {
    
    // 首先调整配置以适应初始条件
    CfdConfig adjusted_config = config;
    adjust_config_for_ic(adjusted_config);
    
    // 使用配置中的参数创建初始条件
    if (!adjusted_config.ic_params.empty()) {
        return create_initial_condition(adjusted_config.ic_type, 
                                       adjusted_config.ic_params);
    } else {
        return create_initial_condition(adjusted_config.ic_type);
    }
}


// ===================================================================
// 工厂函数
// ===================================================================

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

std::unique_ptr<InitialCondition> InitialConditionFactory::create_complex_wave(
    Real a, Real z, Real delta, Real alpha, const std::string& name) {
    return std::make_unique<ComplexWaveInitialCondition>(a, z, delta, alpha, name);
}

std::unique_ptr<InitialCondition> InitialConditionFactory::create_custom(
    const std::function<Real(Real)>& func, const std::string& name) {
    return std::make_unique<CustomInitialCondition>(func, name);
}

// ===================================================================
// 域映射辅助函数
// ===================================================================

bool InitialConditionFactory::requires_specific_domain(const std::string& ic_type) {
    std::string type_lower = to_lower(ic_type);
    
    // 复杂波需要特定域 [-1.0, 1.0]
    return (type_lower == "complex" || type_lower == "complexwave" || 
            type_lower == "precise");
}

std::pair<Real, Real> InitialConditionFactory::get_suggested_domain(const std::string& ic_type) {
    std::string type_lower = to_lower(ic_type);
    
    if (type_lower == "complex" || type_lower == "complexwave" || type_lower == "precise") {
        // 复杂波需要 [-1.0, 1.0] 域以包含所有特征
        return std::make_pair(-1.0, 1.0);
    } else {
        // 默认域 [0.0, 2.0] 适用于大多数测试
        return std::make_pair(0.0, 2.0);
    }
}

void InitialConditionFactory::adjust_config_for_ic(CfdConfig& config) {
    std::string type_lower = to_lower(config.ic_type);
    
    // 检查是否需要调整域
    if (requires_specific_domain(type_lower)) {
        auto suggested_domain = get_suggested_domain(type_lower);
        
        // 如果当前域与建议域不同，打印警告并调整
        if (config.xmin != suggested_domain.first || config.xmax != suggested_domain.second) {
            if (config.verbose) {
                std::cout << "Adjusting domain for " << config.ic_type 
                          << " initial condition:" << std::endl;
                std::cout << "  From: [" << config.xmin << ", " << config.xmax << "]" << std::endl;
                std::cout << "  To:   [" << suggested_domain.first << ", " 
                          << suggested_domain.second << "]" << std::endl;
            }
            
            config.xmin = suggested_domain.first;
            config.xmax = suggested_domain.second;
        }
    }
    
    // 为复杂波设置默认参数（如果未提供）
    if ((type_lower == "complex" || type_lower == "complexwave" || type_lower == "precise") &&
        config.ic_params.empty()) {
        config.ic_params = {0.5, -0.7, 0.005, 10.0};
        if (config.verbose) {
            std::cout << "Using default complex wave parameters: "
                      << "a=0.5, z=-0.7, δ=0.005, α=10.0" << std::endl;
        }
    }
    
    // 为步函数设置默认参数（如果未提供）
    if (type_lower == "step" && config.ic_params.empty()) {
        config.ic_params = {1.0, 2.0, 0.5, 1.0};
    }
    
    // 为正弦波设置默认参数（如果未提供）
    if ((type_lower == "sine" || type_lower == "sinusoidal" || type_lower == "sin") &&
        config.ic_params.empty()) {
        config.ic_params = {0.5, 1.0, 1.0, 0.0};
    }
}

// ===================================================================
// 其他函数
// ===================================================================

std::vector<std::string> InitialConditionFactory::available_initial_conditions() {
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
        
        // 添加域要求信息
        if (requires_specific_domain(ic_type)) {
            auto domain = get_suggested_domain(ic_type);
            oss << "  Suggested domain: [" << domain.first << ", " << domain.second << "]\n";
        }
        
        return oss.str();
    } catch (const std::exception& e) {
        return "Error: " + std::string(e.what());
    }
}

bool InitialConditionFactory::is_available(const std::string& ic_type) {
    std::vector<std::string> available = {
        "step", "sinusoidal", "gaussian", "complex", "complexwave", "precise", "custom",
        "stepfunction", "sin", "sine", "gauss"
    };
    
    std::string type_lower = to_lower(ic_type);
    
    if (std::find(available.begin(), available.end(), type_lower) != available.end()) {
        return true;
    }
    
    // Check custom registry
    auto& registry = get_custom_registry();
    return registry.find(type_lower) != registry.end();
}

std::unique_ptr<InitialCondition> InitialConditionFactory::create_default() {
    // 默认返回复杂波初始条件（与Fortran一致）
    return create_complex_wave();
}

void InitialConditionFactory::register_custom_type(const std::string& name,
                                                  std::function<std::unique_ptr<InitialCondition>()> creator) {
    std::string name_lower = to_lower(name);
    get_custom_registry()[name_lower] = creator;
}

std::map<std::string, std::function<std::unique_ptr<InitialCondition>()>>& 
InitialConditionFactory::get_custom_registry() {
    static std::map<std::string, std::function<std::unique_ptr<InitialCondition>()>> registry;
    return registry;
}

} // namespace cfd