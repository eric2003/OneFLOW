// ==================== src/time_integrator/ThetaMethod.cpp ====================
#include "time_integrator/ThetaMethod.h"
#include "boundary.hpp"
#include <cmath>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

namespace cfd {

ThetaMethod::ThetaMethod(std::unique_ptr<ResidualCalculator> residual_calculator, 
                         Real theta)
    : TimeIntegrator(std::move(residual_calculator))
    , theta_(theta)
    , matrix_assembled_(false)
    , last_dt_(0.0) {
    
    if (theta_ < 0.0 || theta_ > 1.0) {
        throw std::invalid_argument("Theta must be in [0, 1]");
    }
}

void ThetaMethod::step(Vector& u_with_ghosts,
    Real dt,
    const ComputationalDomain& domain) const {

    const int ncells = domain.mesh().ncells();
    const int ist = domain.ist();
    const Real dx = domain.mesh().dx();
    const Real wave_speed = residual_calculator_->wave_speed();

    // 保存旧解（interior cells）
    Vector u_old(ncells);
    for (int i = 0; i < ncells; ++i) {
        u_old[i] = u_with_ghosts[ist + i];
    }

    // 如果是纯显式 (θ=0)，使用简单前向欧拉更新
    if (theta_ == 0.0) {
        Vector residual(ncells);
        residual_calculator_->compute(u_with_ghosts, residual, domain, dx);

        for (int i = 0; i < ncells; ++i) {
            u_with_ghosts[ist + i] += dt * residual[i];
        }
        return;
    }

    // 隐式方法：需要组装矩阵并求解线性系统

    // 检查是否需要重新组装矩阵（dt 变化时重新计算）
    if (!matrix_assembled_ || dt != last_dt_) {
        assemble_matrix_linear_convection(dt, domain);
        matrix_assembled_ = true;
        last_dt_ = dt;
    }

    // ────────────────────────────────────────────────
    // 计算旧时刻的残差：res_old = L u^n
    // 注意：你的 residual_calculator 返回的是 res = - (flux divergence)，即 res = L u^n（对流方程 du/dt = L u）
    Vector residual_old(ncells);
    residual_calculator_->compute(u_with_ghosts, residual_old, domain, dx);

    // 构建右端项：rhs = u^n + (1 - θ) Δt * res_old
    // 这就是 Crank-Nicolson 的正确右端项形式
    Vector rhs(ncells);
    for (int i = 0; i < ncells; ++i) {
        rhs[i] = u_old[i] + (1.0 - theta_) * dt * residual_old[i];
    }

    // 转换为 Eigen 向量
    Eigen::VectorXd rhs_eigen = Eigen::VectorXd::Map(rhs.data(), ncells);

    // ────────────────────────────────────────────────
    // 调试打印：检查 rhs 是否真的更新了
    double max_rhs_diff = 0.0;
    for (int i = 0; i < ncells; ++i) {
        double diff = std::abs(rhs[i] - u_old[i]);
        if (diff > max_rhs_diff) max_rhs_diff = diff;
    }
    std::cout << "Max |rhs - u_old| = " << max_rhs_diff << "\n";

    // 打印几个典型点，便于观察
    if (ncells > 0) {
        std::cout << "  rhs[0] = " << rhs[0] << ", u_old[0] = " << u_old[0] << "\n";
        std::cout << "  rhs[" << ncells/2 << "] = " << rhs[ncells/2] 
            << ", u_old[" << ncells/2 << "] = " << u_old[ncells/2] << "\n";
    }

    // 应用边界条件到右端项（如果需要）
    apply_boundary_to_system(rhs, domain);

    // ────────────────────────────────────────────────
    // 求解线性系统：(I - θ Δt L) u_new = rhs
    Eigen::BiCGSTAB<Eigen::SparseMatrix<Real>> solver;
    solver.compute(system_matrix_);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Matrix decomposition failed");
    }

    Eigen::VectorXd u_new_vec = solver.solve(rhs_eigen);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Linear solve failed");
    }

    // 打印求解器收敛信息（强烈建议保留）
    std::cout << "BiCGSTAB iterations: " << solver.iterations() 
        << ", estimated error: " << solver.error() << "\n";

    // 更新解（只更新内部格子）
    for (int i = 0; i < ncells; ++i) {
        u_with_ghosts[ist + i] = u_new_vec[i];
    }

    // 应用周期边界条件更新鬼格
    domain.apply_periodic_boundary(u_with_ghosts);
}

//void ThetaMethod::step(Vector& u_with_ghosts, 
//                      Real dt,
//                      const ComputationalDomain& domain) const {
//    
//    const int ncells = domain.mesh().ncells();
//    const int ist = domain.ist();
//    const Real dx = domain.mesh().dx();
//    const Real wave_speed = residual_calculator_->wave_speed();
//    
//    // 保存旧解
//    Vector u_old(ncells);
//    for (int i = 0; i < ncells; ++i) {
//        u_old[i] = u_with_ghosts[ist + i];
//    }
//    
//    // 如果是纯显式 (θ=0)，使用简单更新
//    if (theta_ == 0.0) {
//        // 显式欧拉
//        Vector residual(ncells);
//        residual_calculator_->compute(u_with_ghosts, residual, domain, dx);
//        
//        for (int i = 0; i < ncells; ++i) {
//            u_with_ghosts[ist + i] += dt * residual[i];
//        }
//        
//        return;
//    }
//    
//    // 隐式方法：需要组装和求解线性系统
//    
//    // 检查是否需要重新组装矩阵
//    if (!matrix_assembled_ || dt != last_dt_) {
//        assemble_matrix_linear_convection(dt, domain);
//        matrix_assembled_ = true;
//        last_dt_ = dt;
//    }
//    
//    // 计算右端项: rhs = [I + (1-θ)dt * L] * u_old
//    Eigen::VectorXd u_old_vec = Eigen::VectorXd::Map(u_old.data(), ncells);
//    Eigen::VectorXd rhs_vec = system_matrix_ * u_old_vec;  // 注意：这里system_matrix_是[I - θdtL]
//    
//    // 转换为Vector
//    Vector rhs(ncells);
//    for (int i = 0; i < ncells; ++i) {
//        rhs[i] = rhs_vec[i];
//    }
//    
//    // 应用边界条件到右端项
//    apply_boundary_to_system(rhs, domain);
//    
//    // 求解线性系统: [I - θdtL] * u_new = rhs
//    Eigen::BiCGSTAB<Eigen::SparseMatrix<Real>> solver;
//    solver.compute(system_matrix_);
//    
//    if (solver.info() != Eigen::Success) {
//        throw std::runtime_error("Matrix decomposition failed");
//    }
//
//    // ★★★ 这里就是 “solve 之前” ★★★
//    // 在调用 solver.solve() 之前，rhs_eigen 已经准备好了
//    // 你可以在这里插入打印 rhs 和 u_old 的差
//
//    // 打印 max |rhs - u_old| 的代码就加在这里
//    double max_rhs_diff = 0.0;
//    for(int i = 0; i < ncells; ++i) {
//        double diff = std::abs(rhs[i] - u_old[i]);  // rhs 是你的 std::vector
//        if (diff > max_rhs_diff) max_rhs_diff = diff;
//    }
//    std::cout << "Max |rhs - u_old| = " << max_rhs_diff << "\n";
//
//    // 也可以打印几个典型位置的值，便于观察
//    std::cout << "rhs[0] = " << rhs[0] << ", u_old[0] = " << u_old[0] << "\n";
//    std::cout << "rhs[100] = " << rhs[100] << ", u_old[100] = " << u_old[100] << "\n";
//    
//    Eigen::VectorXd rhs_eigen = Eigen::VectorXd::Map(rhs.data(), ncells);
//    Eigen::VectorXd u_new_vec = solver.solve(rhs_eigen);
//    
//    if (solver.info() != Eigen::Success) {
//        throw std::runtime_error("Linear solve failed");
//    }
//    
//    // 更新解
//    for (int i = 0; i < ncells; ++i) {
//        u_with_ghosts[ist + i] = u_new_vec[i];
//    }
//    
//    // 应用边界条件更新虚网格
//    domain.apply_periodic_boundary(u_with_ghosts);
//}

//void ThetaMethod::assemble_matrix_linear_convection(Real dt, 
//                                                   const ComputationalDomain& domain) const {
//    
//    const int ncells = domain.mesh().ncells();
//    const Real dx = domain.mesh().dx();
//    const Real wave_speed = residual_calculator_->wave_speed();
//    const Real a = wave_speed;
//    
//    // 使用三对角矩阵（针对中心差分空间离散）
//    // 对于线性对流方程 du/dt + a * du/dx = 0
//    // 中心差分: du/dx ≈ (u_{i+1} - u_{i-1}) / (2dx)
//    // 所以 L * u = -a * (u_{i+1} - u_{i-1}) / (2dx)
//    
//    std::vector<Eigen::Triplet<Real>> triplets;
//    triplets.reserve(3 * ncells);
//    
//    const Real alpha = -theta_ * dt * a / (2.0 * dx);
//    
//    for (int i = 0; i < ncells; ++i) {
//        // 主对角线: 1 - θdt * L_ii
//        triplets.push_back(Eigen::Triplet<Real>(i, i, 1.0));
//        
//        // 下对角线: -θdt * L_{i,i-1}
//        if (i > 0) {
//            triplets.push_back(Eigen::Triplet<Real>(i, i-1, -alpha));
//        }
//        
//        // 上对角线: -θdt * L_{i,i+1}
//        if (i < ncells - 1) {
//            triplets.push_back(Eigen::Triplet<Real>(i, i+1, alpha));
//        }
//    }
//    
//    // 构建稀疏矩阵
//    system_matrix_.resize(ncells, ncells);
//    system_matrix_.setFromTriplets(triplets.begin(), triplets.end());
//    system_matrix_.makeCompressed();
//
//    // 在 assemble_matrix_linear_convection 最后加
//    std::cout << "\n=== System Matrix (first few rows) ===\n";
//    for (int i = 0; i < std::min(5, ncells); ++i) {
//        std::cout << "Row " << i << ": ";
//        for (Eigen::SparseMatrix<double>::InnerIterator it(system_matrix_, i); it; ++it) {
//            std::cout << "  col=" << it.col() << " val=" << it.value();
//        }
//        std::cout << "\n";
//    }
//}


void ThetaMethod::assemble_matrix_linear_convection(Real dt,
    const ComputationalDomain& domain) const {
    const int n = domain.mesh().ncells();
    const Real a = residual_calculator_->wave_speed();  // 注意符号
    const Real dx = domain.mesh().dx();

    // alpha = θ * dt * a / (2 dx)
    // 中心差分：f_{i+1/2} = 0.5 a (u_i + u_{i+1})
    // → L u = -a (u_{i+1} - u_{i-1}) / (2 dx)
    const Real alpha = theta_ * dt * a / (2.0 * dx);


    std::cout << "theta=" << theta_ << ", dt=" << dt 
        << ", a=" << a << ", dx=" << dx 
        << " → alpha = " << alpha << "\n";

    std::vector<Eigen::Triplet<Real>> triplets;
    triplets.reserve(3 * n);  // 每行最多 3 个非零

    for (int i = 0; i < n; ++i) {
        // 主对角：永远只加一次 1.0
        triplets.emplace_back(i, i, 1.0);

        // 左邻居贡献：- θ dt a /(2 dx)
        int left = (i == 0) ? n-1 : i-1;           // 周期
        triplets.emplace_back(i, left, -alpha);

        // 右邻居贡献：+ θ dt a /(2 dx)
        int right = (i == n-1) ? 0 : i+1;          // 周期
        triplets.emplace_back(i, right, +alpha);
    }

    system_matrix_.resize(n, n);
    system_matrix_.setFromTriplets(triplets.begin(), triplets.end());

    std::cout << "\nDiagonal elements (first 5):\n";
    for(int i=0; i<5; i++){
        std::cout << "A(" << i << "," << i << ") = " << system_matrix_.coeff(i,i) << "\n";
    }
    system_matrix_.makeCompressed();

    // 调试打印（建议保留一段时间）
    std::cout << "Matrix assembled. nnz = " << system_matrix_.nonZeros() << "\n";
}

std::string ThetaMethod::name() const {
    if (theta_ == 0.0) {
        return "Explicit Euler (θ=0.0)";
    } else if (theta_ == 0.5) {
        return "Crank-Nicolson (θ=0.5)";
    } else if (theta_ == 1.0) {
        return "Implicit Euler (θ=1.0)";
    } else {
        return "Theta Method (θ=" + std::to_string(theta_) + ")";
    }
}

int ThetaMethod::order() const {
    if (theta_ == 0.5) {
        return 2;  // Crank-Nicolson是二阶
    } else {
        return 1;  // 欧拉方法是一阶
    }
}

void ThetaMethod::apply_boundary_to_system(Vector& rhs,
                                          const ComputationalDomain& domain) const {
    // 对于周期边界条件，系统矩阵已经是循环三对角的
    // 不需要特殊处理右端项
    // 如果是其他边界条件，需要在这里调整
}

} // namespace cfd