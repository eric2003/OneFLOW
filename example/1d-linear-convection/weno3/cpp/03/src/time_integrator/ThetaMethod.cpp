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
    
    // 保存旧解
    Vector u_old(ncells);
    for (int i = 0; i < ncells; ++i) {
        u_old[i] = u_with_ghosts[ist + i];
    }
    
    // 如果是纯显式 (θ=0)，使用简单更新
    if (theta_ == 0.0) {
        // 显式欧拉
        Vector residual(ncells);
        residual_calculator_->compute(u_with_ghosts, residual, domain, dx);
        
        for (int i = 0; i < ncells; ++i) {
            u_with_ghosts[ist + i] += dt * residual[i];
        }
        
        return;
    }
    
    // 隐式方法：需要组装和求解线性系统
    
    // 检查是否需要重新组装矩阵
    if (!matrix_assembled_ || dt != last_dt_) {
        assemble_matrix_linear_convection(dt, domain);
        matrix_assembled_ = true;
        last_dt_ = dt;
    }
    
    // 计算右端项: rhs = [I + (1-θ)dt * L] * u_old
    Eigen::VectorXd u_old_vec = Eigen::VectorXd::Map(u_old.data(), ncells);
    Eigen::VectorXd rhs_vec = system_matrix_ * u_old_vec;  // 注意：这里system_matrix_是[I - θdtL]
    
    // 转换为Vector
    Vector rhs(ncells);
    for (int i = 0; i < ncells; ++i) {
        rhs[i] = rhs_vec[i];
    }
    
    // 应用边界条件到右端项
    apply_boundary_to_system(rhs, domain);
    
    // 求解线性系统: [I - θdtL] * u_new = rhs
    Eigen::BiCGSTAB<Eigen::SparseMatrix<Real>> solver;
    solver.compute(system_matrix_);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Matrix decomposition failed");
    }
    
    Eigen::VectorXd rhs_eigen = Eigen::VectorXd::Map(rhs.data(), ncells);
    Eigen::VectorXd u_new_vec = solver.solve(rhs_eigen);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Linear solve failed");
    }
    
    // 更新解
    for (int i = 0; i < ncells; ++i) {
        u_with_ghosts[ist + i] = u_new_vec[i];
    }
    
    // 应用边界条件更新虚网格
    domain.apply_periodic_boundary(u_with_ghosts);
}

void ThetaMethod::assemble_matrix_linear_convection(Real dt, 
                                                   const ComputationalDomain& domain) const {
    
    const int ncells = domain.mesh().ncells();
    const Real dx = domain.mesh().dx();
    const Real wave_speed = residual_calculator_->wave_speed();
    const Real a = wave_speed;
    
    // 使用三对角矩阵（针对中心差分空间离散）
    // 对于线性对流方程 du/dt + a * du/dx = 0
    // 中心差分: du/dx ≈ (u_{i+1} - u_{i-1}) / (2dx)
    // 所以 L * u = -a * (u_{i+1} - u_{i-1}) / (2dx)
    
    std::vector<Eigen::Triplet<Real>> triplets;
    triplets.reserve(3 * ncells);
    
    const Real alpha = -theta_ * dt * a / (2.0 * dx);
    
    for (int i = 0; i < ncells; ++i) {
        // 主对角线: 1 - θdt * L_ii
        triplets.push_back(Eigen::Triplet<Real>(i, i, 1.0));
        
        // 下对角线: -θdt * L_{i,i-1}
        if (i > 0) {
            triplets.push_back(Eigen::Triplet<Real>(i, i-1, -alpha));
        }
        
        // 上对角线: -θdt * L_{i,i+1}
        if (i < ncells - 1) {
            triplets.push_back(Eigen::Triplet<Real>(i, i+1, alpha));
        }
    }
    
    // 构建稀疏矩阵
    system_matrix_.resize(ncells, ncells);
    system_matrix_.setFromTriplets(triplets.begin(), triplets.end());
    system_matrix_.makeCompressed();
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