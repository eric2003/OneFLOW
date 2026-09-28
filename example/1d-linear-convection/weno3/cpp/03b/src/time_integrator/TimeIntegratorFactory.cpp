// ==================== src/time_integrator/TimeIntegratorFactory.cpp ====================
#include "time_integrator/TimeIntegratorFactory.h"
#include "time_integrator/ThetaMethod.h"  // 新增

namespace cfd {

// 主创建函数（处理所有方法）
std::unique_ptr<TimeIntegrator> TimeIntegratorFactory::create_integrator(
    const std::string& method,
    std::unique_ptr<ResidualCalculator> residual_calculator) {
    
    std::string method_lower = method;
    std::transform(method_lower.begin(), method_lower.end(), 
                   method_lower.begin(), ::tolower);
    
    // 检查方法类型并分发
    if (method_lower.find("rk") == 0 || 
        method_lower == "euler" || 
        method_lower == "forwardeuler" ||
        method_lower == "midpoint") {
        // 显式方法
        return create_explicit_integrator(method_lower, std::move(residual_calculator));
    } 
    else if (method_lower.find("crank") != std::string::npos ||
             method_lower.find("implicit") != std::string::npos ||
             method_lower.find("theta") != std::string::npos) {
        // 隐式方法
        return create_implicit_integrator(method_lower, std::move(residual_calculator));
    }
    else {
        throw std::invalid_argument("Unknown time integration method: " + method);
    }
}

// 创建显式RK积分器
std::unique_ptr<TimeIntegrator> TimeIntegratorFactory::create_explicit_integrator(
    const std::string& method,
    std::unique_ptr<ResidualCalculator> residual_calculator) {
    
    std::string method_lower = method;
    std::transform(method_lower.begin(), method_lower.end(), 
                   method_lower.begin(), ::tolower);
    
    if (method_lower == "rk1" || method_lower == "euler" || 
        method_lower == "forwardeuler") {
        return std::make_unique<RK1Integrator>(std::move(residual_calculator));
    }
    else if (method_lower == "rk2" || method_lower == "midpoint") {
        return std::make_unique<RK2Integrator>(std::move(residual_calculator));
    }
    else if (method_lower == "rk3" || method_lower == "tvdrk3") {
        return std::make_unique<RK3Integrator>(std::move(residual_calculator));
    }
    else {
        // 尝试按阶数解析
        try {
            int order = std::stoi(method_lower.substr(2)); // 提取"rk"后面的数字
            return create_explicit_integrator(order, std::move(residual_calculator));
        } catch (...) {
            throw std::invalid_argument("Unknown explicit method: " + method);
        }
    }
}

// 创建显式RK积分器（按阶数） - 保持向后兼容
std::unique_ptr<TimeIntegrator> TimeIntegratorFactory::create_integrator(
    int order,
    std::unique_ptr<ResidualCalculator> residual_calculator) {
    
    // 这个函数现在明确只处理显式RK方法
    return create_explicit_integrator(order, std::move(residual_calculator));
}

// 显式RK积分器（按阶数）的内部实现
std::unique_ptr<TimeIntegrator> TimeIntegratorFactory::create_explicit_integrator(
    int order,
    std::unique_ptr<ResidualCalculator> residual_calculator) {
    
    switch (order) {
        case 1:
            return std::make_unique<RK1Integrator>(std::move(residual_calculator));
        case 2:
            return std::make_unique<RK2Integrator>(std::move(residual_calculator));
        case 3:
            return std::make_unique<RK3Integrator>(std::move(residual_calculator));
        default:
            throw std::invalid_argument("Unsupported RK order: " + std::to_string(order) +
                                      " (only 1, 2, 3 are supported)");
    }
}
// 创建隐式积分器
std::unique_ptr<TimeIntegrator> TimeIntegratorFactory::create_implicit_integrator(
    const std::string& method,
    std::unique_ptr<ResidualCalculator> residual_calculator,
    Real theta) {
    
    std::string method_lower = method;
    std::transform(method_lower.begin(), method_lower.end(), 
                   method_lower.begin(), ::tolower);
    
    if (method_lower == "crank-nicolson" || method_lower == "cn") {
        return std::make_unique<ThetaMethod>(std::move(residual_calculator), 0.5);
    }
    else if (method_lower == "implicit-euler" || method_lower == "backward-euler" || 
             method_lower == "be") {
        return std::make_unique<ThetaMethod>(std::move(residual_calculator), 1.0);
    }
    else if (method_lower == "explicit-euler") {
        return std::make_unique<ThetaMethod>(std::move(residual_calculator), 0.0);
    }
    else if (method_lower.find("theta") != std::string::npos) {
        // 解析theta值
        Real actual_theta = theta; // 使用默认或传入值
        
        if (method_lower.find("theta-") == 0) {
            try {
                std::string theta_str = method_lower.substr(6);
                actual_theta = std::stod(theta_str);
            } catch (...) {
                // 保持默认值
            }
        }
        
        // 验证theta范围
        if (actual_theta < 0.0 || actual_theta > 1.0) {
            throw std::invalid_argument("Theta must be in [0, 1], got: " + 
                                      std::to_string(actual_theta));
        }
        
        return std::make_unique<ThetaMethod>(std::move(residual_calculator), actual_theta);
    }
    else {
        throw std::invalid_argument("Unknown implicit method: " + method);
    }
}

// 检查方法是否可用
bool TimeIntegratorFactory::is_available(const std::string& method) {
    std::string method_lower = method;
    std::transform(method_lower.begin(), method_lower.end(), 
                   method_lower.begin(), ::tolower);
    
    // 显式方法
    if (method_lower == "rk1" || method_lower == "euler" || 
        method_lower == "forwardeuler" ||
        method_lower == "rk2" || method_lower == "midpoint" ||
        method_lower == "rk3" || method_lower == "tvdrk3") {
        return true;
    }
    
    // 隐式方法
    if (method_lower == "crank-nicolson" || method_lower == "cn" ||
        method_lower == "implicit-euler" || method_lower == "backward-euler" ||
        method_lower == "be" || method_lower == "explicit-euler" ||
        method_lower == "theta" || method_lower == "theta-0.25" ||
        method_lower == "theta-0.5" || method_lower == "theta-0.75" ||
        method_lower == "theta-1.0") {
        return true;
    }
    
    // 检查是否theta-<value>格式
    if (method_lower.find("theta-") == 0) {
        try {
            std::string theta_str = method_lower.substr(6);
            Real theta = std::stod(theta_str);
            return (theta >= 0.0 && theta <= 1.0);
        } catch (...) {
            return false;
        }
    }
    
    // 检查是否rk<数字>格式
    if (method_lower.find("rk") == 0) {
        try {
            int order = std::stoi(method_lower.substr(2));
            return (order >= 1 && order <= 3);
        } catch (...) {
            return false;
        }
    }
    
    return false;
}

// 获取所有可用方法
std::vector<std::string> TimeIntegratorFactory::available_integrators() {
    return {
        "rk1 (Forward Euler) - 1st order explicit",
        "rk2 (Midpoint) - 2nd order explicit", 
        "rk3 (TVD RK3) - 3rd order explicit",
        "crank-nicolson - 2nd order implicit (θ=0.5)",
        "implicit-euler - 1st order implicit (θ=1.0)",
        "explicit-euler - 1st order explicit (θ=0.0)",
        "theta-0.25 - Theta method with θ=0.25",
        "theta-0.75 - Theta method with θ=0.75",
        "theta - General theta method (specify θ value)"
    };
}

// 按类型获取方法
std::vector<std::string> TimeIntegratorFactory::available_explicit_integrators() {
    return {"rk1", "rk2", "rk3", "euler", "forwardeuler", "midpoint", "tvdrk3"};
}

std::vector<std::string> TimeIntegratorFactory::available_implicit_integrators() {
    return {"crank-nicolson", "cn", "implicit-euler", "backward-euler", 
            "be", "explicit-euler", "theta", "theta-0.25", "theta-0.5", 
            "theta-0.75", "theta-1.0"};
}


std::string TimeIntegratorFactory::get_integrator_info(const std::string& method) {
    // 新增θ方法信息
    if (method == "crank-nicolson" || method == "CN") {
        return "Crank-Nicolson: 2nd order implicit, unconditionally stable";
    }
    else if (method == "implicit-euler" || method == "backward-euler") {
        return "Implicit Euler: 1st order implicit, unconditionally stable";
    }
    else if (method.find("theta") != std::string::npos) {
        return "Theta method: general implicit-explicit, θ in [0,1]";
    }
    else  if (method == "rk1" || method == "euler" || method == "ForwardEuler" || method == "RK1") {
        return "Forward Euler (RK1): 1st order explicit, CFL limit = 1.0";
    }
    else if (method == "rk2" || method == "RK2" || method == "midpoint") {
        return "Runge-Kutta 2 (RK2): 2nd order explicit, CFL limit = 1.0";
    }
    else if (method == "rk3" || method == "RK3" || method == "TVDRK3") {
        return "TVD Runge-Kutta 3 (RK3): 3rd order explicit, TVD preserving, CFL limit = 1.0";
    }
    else {
        return "Unknown time integration method: " + method;
    }
}

std::unique_ptr<TimeIntegrator> TimeIntegratorFactory::create_default_integrator(
    std::unique_ptr<ResidualCalculator> residual_calculator) {
    return std::make_unique<RK2Integrator>(std::move(residual_calculator));
}

} // namespace cfd