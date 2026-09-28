#include "time_integrator/TimeIntegratorFactory.h"

namespace cfd {

std::unique_ptr<TimeIntegrator> TimeIntegratorFactory::create_integrator(
    const std::string& method,
    std::unique_ptr<ResidualCalculator> residual_calculator) {
    
    if (method == "rk1" || method == "euler" || method == "ForwardEuler" || method == "RK1") {
        return std::make_unique<RK1Integrator>(std::move(residual_calculator));
    }
    else if (method == "rk2" || method == "RK2" || method == "midpoint") {
        return std::make_unique<RK2Integrator>(std::move(residual_calculator));
    }
    else if (method == "rk3" || method == "RK3" || method == "TVDRK3") {
        return std::make_unique<RK3Integrator>(std::move(residual_calculator));
    }
    else {
        throw std::invalid_argument("Unknown time integration method: " + method);
    }
}

std::unique_ptr<TimeIntegrator> TimeIntegratorFactory::create_integrator(
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
            throw std::invalid_argument("Unsupported RK order: " + std::to_string(order));
    }
}

std::unique_ptr<TimeIntegrator> TimeIntegratorFactory::create_integrator(
    const CfdConfig& config,
    std::unique_ptr<ResidualCalculator> residual_calculator) {
    
    return create_integrator(config.rk_order, std::move(residual_calculator));
}

std::vector<std::string> TimeIntegratorFactory::available_integrators() {
    return {
        "rk1 - Forward Euler (1st order)",
        "rk2 - Runge-Kutta 2 (2nd order)",
        "rk3 - TVD Runge-Kutta 3 (3rd order)"
    };
}

std::string TimeIntegratorFactory::get_integrator_info(const std::string& method) {
    if (method == "rk1" || method == "euler" || method == "ForwardEuler" || method == "RK1") {
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

bool TimeIntegratorFactory::is_available(const std::string& method) {
    std::vector<std::string> available = {
        "rk1", "euler", "ForwardEuler", "RK1",
        "rk2", "RK2", "midpoint",
        "rk3", "RK3", "TVDRK3"
    };
    
    return std::find(available.begin(), available.end(), method) != available.end();
}

std::unique_ptr<TimeIntegrator> TimeIntegratorFactory::create_default_integrator(
    std::unique_ptr<ResidualCalculator> residual_calculator) {
    return std::make_unique<RK2Integrator>(std::move(residual_calculator));
}

} // namespace cfd