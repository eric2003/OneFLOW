#ifndef INITIAL_CONDITION_FACTORY_H
#define INITIAL_CONDITION_FACTORY_H

#include "InitialCondition.h"
#include "StepInitialCondition.h"
#include "SinusoidalInitialCondition.h"
#include "GaussianInitialCondition.h"
#include "CustomInitialCondition.h"
#include "ComplexWaveInitialCondition.h"  // 添加复杂波头文件
#include "config.hpp"
#include <memory>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

namespace cfd {

// ===================================================================
// Initial condition factory
// ===================================================================
class InitialConditionFactory {
public:
    // Create initial condition by type name with parameters
    static std::unique_ptr<InitialCondition> create_initial_condition(
        const std::string& ic_type,
        const std::vector<Real>& params = {});
    
    // Create initial condition from configuration - 新增方法
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
    
    // Create complex wave (from Fortran)
    static std::unique_ptr<InitialCondition> create_complex_wave(
        Real a = 0.5, Real z = -0.7, Real delta = 0.005,
        Real alpha = 10.0, const std::string& name = "Complex Wave");
    
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
    
    // Create default initial condition
    static std::unique_ptr<InitialCondition> create_default();
    
    // Register custom initial condition type
    static void register_custom_type(const std::string& name,
                                   std::function<std::unique_ptr<InitialCondition>()> creator);
    
    // ===================================================================
    // 域映射辅助函数（用于复杂波等需要特定域的初始条件）
    // ===================================================================
    
    // 检查初始条件是否需要特定域
    static bool requires_specific_domain(const std::string& ic_type);
    
    // 获取建议的域范围
    static std::pair<Real, Real> get_suggested_domain(const std::string& ic_type);
    
    // 调整配置以适应初始条件
    static void adjust_config_for_ic(CfdConfig& config);
    
private:
    // Registry for custom initial conditions
    static std::map<std::string, std::function<std::unique_ptr<InitialCondition>()>>& 
    get_custom_registry();
    
    // 将字符串转换为小写
    static std::string to_lower(const std::string& str) {
        std::string result = str;
        std::transform(result.begin(), result.end(), result.begin(), ::tolower);
        return result;
    }
};

} // namespace cfd

#endif // INITIAL_CONDITION_FACTORY_H