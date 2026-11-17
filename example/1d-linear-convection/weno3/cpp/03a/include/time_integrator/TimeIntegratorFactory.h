// ==================== include/time_integrator/TimeIntegratorFactory.h ====================
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

    class TimeIntegratorFactory {
    public:
        // ===================================================================
        // 主创建函数
        // ===================================================================

        // 通过字符串创建（支持所有方法）
        static std::unique_ptr<TimeIntegrator> create_integrator(
            const std::string& method,
            std::unique_ptr<ResidualCalculator> residual_calculator);

        // 通过阶数创建（保持向后兼容，只处理显式RK方法）
        static std::unique_ptr<TimeIntegrator> create_integrator(
            int order,
            std::unique_ptr<ResidualCalculator> residual_calculator);

        // ===================================================================
        // 显式方法创建函数
        // ===================================================================

        // 通过字符串创建显式方法
        static std::unique_ptr<TimeIntegrator> create_explicit_integrator(
            const std::string& method,
            std::unique_ptr<ResidualCalculator> residual_calculator);

        // 通过阶数创建显式RK方法
        static std::unique_ptr<TimeIntegrator> create_explicit_integrator(
            int order,
            std::unique_ptr<ResidualCalculator> residual_calculator);

        // ===================================================================
        // 隐式方法创建函数
        // ===================================================================

        // 通过字符串创建隐式方法
        static std::unique_ptr<TimeIntegrator> create_implicit_integrator(
            const std::string& method,
            std::unique_ptr<ResidualCalculator> residual_calculator,
            Real theta = 0.5);

        // ===================================================================
        // 其他功能函数
        // ===================================================================

        // 从配置创建
        static std::unique_ptr<TimeIntegrator> create_integrator(
            const CfdConfig& config,
            std::unique_ptr<ResidualCalculator> residual_calculator);

        // 检查方法是否可用
        static bool is_available(const std::string& method);

        // 获取所有可用方法
        static std::vector<std::string> available_integrators();

        // 按类型获取方法
        static std::vector<std::string> available_explicit_integrators();
        static std::vector<std::string> available_implicit_integrators();

        // 获取方法信息
        static std::string get_integrator_info(const std::string& method);

        // 获取默认积分器
        static std::unique_ptr<TimeIntegrator> create_default_integrator(
            std::unique_ptr<ResidualCalculator> residual_calculator);
    };

} // namespace cfd

#endif // TIME_INTEGRATOR_FACTORY_H
