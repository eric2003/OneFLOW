// ==================== include/analysis.hpp (简化版) ====================
#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include "cfd_solver.hpp"
#include <vector>
#include <string>
#include <memory>

namespace cfd {

    // ===================================================================
    // 分析工具类（纯工具函数，不包含具体配置）
    // ===================================================================
    class Analysis {
    public:
        // ===================================================================
        // 通用分析工具（接受配置作为参数）
        // ===================================================================

        // 运行两个方案的比较
        static void compare_schemes(const CfdConfig& config1,
            const CfdConfig& config2,
            bool verbose = true);

        // 运行收敛性研究
        static void run_convergence_study(const CfdConfig& base_config,
            const std::vector<int>& cell_counts);

        // 运行通量比较
        static void compare_fluxes(const CfdConfig& base_config,
            const std::vector<std::string>& flux_types);

        // 运行时间积分器比较
        static void compare_time_integrators(const CfdConfig& base_config,
            const std::vector<std::string>& methods);

        // ===================================================================
        // 工具函数（纯数学计算，最稳定）
        // ===================================================================
        static double compute_l1_error(const Vector& numerical, const Vector& analytical);
        static double compute_l2_error(const Vector& numerical, const Vector& analytical);
        static double compute_linf_error(const Vector& numerical, const Vector& analytical);

        // ===================================================================
        // 结果处理工具
        // ===================================================================
        static void write_comparison_results(const std::string& filename,
            const Vector& x,
            const Vector& u1,
            const Vector& u2,
            const Vector& analytical);

        static void plot_results(const std::string& output_dir,
            const std::vector<std::string>& labels,
            const std::vector<Vector>& solutions);

    private:
        // 禁止实例化（纯工具类）
        Analysis() = delete;
        ~Analysis() = delete;
    };

} // namespace cfd

#endif // ANALYSIS_HPP