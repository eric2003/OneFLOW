// ==================== src/core/analysis.cpp (简化版) ====================
#include "analysis.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <sstream>

namespace cfd {

    // 纯工具函数实现
    double Analysis::compute_l1_error(const Vector& numerical, const Vector& analytical) {
        if (numerical.size() != analytical.size()) {
            throw std::invalid_argument("Vector sizes do not match");
        }

        double error = 0.0;
        for (size_t i = 0; i < numerical.size(); ++i) {
            error += std::abs(numerical[i] - analytical[i]);
        }
        return error / numerical.size();
    }

    double Analysis::compute_l2_error(const Vector& numerical, const Vector& analytical) {
        if (numerical.size() != analytical.size()) {
            throw std::invalid_argument("Vector sizes do not match");
        }

        double error = 0.0;
        for (size_t i = 0; i < numerical.size(); ++i) {
            double diff = numerical[i] - analytical[i];
            error += diff * diff;
        }
        return std::sqrt(error / numerical.size());
    }

    double Analysis::compute_linf_error(const Vector& numerical, const Vector& analytical) {
        if (numerical.size() != analytical.size()) {
            throw std::invalid_argument("Vector sizes do not match");
        }

        double max_error = 0.0;
        for (size_t i = 0; i < numerical.size(); ++i) {
            double error = std::abs(numerical[i] - analytical[i]);
            if (error > max_error) max_error = error;
        }
        return max_error;
    }

    // 通用比较函数
    void Analysis::compare_schemes(const CfdConfig& config1,
        const CfdConfig& config2,
        bool verbose) {
        if (verbose) {
            std::cout << "==========================================" << std::endl;
            std::cout << "Scheme Comparison" << std::endl;
            std::cout << "==========================================" << std::endl;
            std::cout << "Scheme 1: " << config1.recon_scheme << config1.spatial_order << std::endl;
            std::cout << "Scheme 2: " << config2.recon_scheme << config2.spatial_order << std::endl;
        }

        // 运行模拟
        Vector u1 = CfdSolver::run_single_simulation(config1);
        Vector u2 = CfdSolver::run_single_simulation(config2);

        // 计算解析解
        CfdSolver solver(config1);
        Vector u_analytical = solver.compute_analytical_solution(config1.final_time);

        // 计算误差
        double error1 = compute_l1_error(u1, u_analytical);
        double error2 = compute_l1_error(u2, u_analytical);

        if (verbose) {
            std::cout << "\nError analysis (L1 norm):" << std::endl;
            std::cout << "  " << config1.recon_scheme << ": " << error1 << std::endl;
            std::cout << "  " << config2.recon_scheme << ": " << error2 << std::endl;
            std::cout << "  Error ratio: " << error1 / error2 << std::endl;
        }

        // 写入结果
        const auto& mesh = solver.domain().mesh();
        write_comparison_results("scheme_comparison.txt",
            mesh.cell_centers(),
            u1, u2, u_analytical);
    }

    void Analysis::run_convergence_study(const CfdConfig& base_config,
        const std::vector<int>& cell_counts) {
        std::cout << "==========================================" << std::endl;
        std::cout << "Convergence Study" << std::endl;
        std::cout << "==========================================" << std::endl;

        std::vector<double> errors;
        std::vector<double> resolutions;

        for (int cells : cell_counts) {
            // 复制配置并修改网格
            CfdConfig config = base_config;
            config.ncells = cells;

            std::cout << "\nRunning with " << cells << " cells..." << std::endl;

            // 运行模拟
            Vector u_numerical = CfdSolver::run_single_simulation(config);

            // 计算解析解
            CfdSolver solver(config);
            Vector u_analytical = solver.compute_analytical_solution(config.final_time);

            // 计算误差
            double error = compute_l1_error(u_numerical, u_analytical);
            errors.push_back(error);
            resolutions.push_back(1.0 / cells);

            std::cout << "  L1 error: " << error << std::endl;
        }

        // 计算收敛率
        std::cout << "\nConvergence rates:" << std::endl;
        for (size_t i = 1; i < errors.size(); ++i) {
            double rate = std::log(errors[i-1] / errors[i]) / 
                std::log(resolutions[i-1] / resolutions[i]);
            std::cout << "  " << cell_counts[i-1] << " → " << cell_counts[i] 
                << ": " << rate << std::endl;
        }
    }

    // 工具函数实现
    void Analysis::write_comparison_results(const std::string& filename,
        const Vector& x,
        const Vector& u1,
        const Vector& u2,
        const Vector& analytical) {
        std::ofstream file(filename);
        file << "# Scheme Comparison Results" << std::endl;
        file << "# X, Solution1, Solution2, Analytical" << std::endl;

        for (size_t i = 0; i < x.size(); ++i) {
            file << std::setw(12) << x[i]
                << std::setw(12) << u1[i]
                << std::setw(12) << u2[i]
                << std::setw(12) << analytical[i] << std::endl;
        }

        file.close();
    }

} // namespace cfd