// ==================== include/time_integrator/ThetaMethod.h ====================
#ifndef THETA_METHOD_H
#define THETA_METHOD_H

#include "TimeIntegrator.h"
#include <Eigen/Sparse>

namespace cfd {

/**
 * @brief θ-方法统一框架
 * 
 * θ = 0.0: 显式欧拉 (explicit Euler)
 * θ = 0.5: Crank-Nicolson
 * θ = 1.0: 隐式欧拉 (implicit Euler)
 */
class ThetaMethod : public TimeIntegrator {
private:
    Real theta_;  // θ ∈ [0, 1]
    
    // 线性求解器状态
    mutable Eigen::SparseMatrix<Real> system_matrix_;
    mutable Eigen::VectorXd rhs_;
    mutable bool matrix_assembled_;
    mutable Real last_dt_;
    
    // 三对角求解器（针对中心差分）
    void solve_tridiagonal(const Vector& a, const Vector& b, 
                          const Vector& c, const Vector& d,
                          Vector& x) const;
    
public:
    ThetaMethod(std::unique_ptr<ResidualCalculator> residual_calculator, 
                Real theta = 0.5);
    ~ThetaMethod() override = default;
    
    void step(Vector& u_with_ghosts, 
              Real dt,
              const ComputationalDomain& domain) const override;
    
    std::string name() const override;
    int order() const override;
    
    bool is_explicit() const override { return (theta_ == 0.0); }
    bool is_implicit() const override { return (theta_ > 0.0); }
    
    // Crank-Nicolson (θ=0.5) 是无条件稳定的
    bool is_unconditionally_stable() const { return (theta_ >= 0.5); }
    
    // θ参数访问器
    Real theta() const { return theta_; }
    
private:
    // 针对线性对流方程的矩阵组装（使用中心差分空间离散）
    void assemble_matrix_linear_convection(Real dt, 
                                          const ComputationalDomain& domain) const;
    
    // 通用矩阵组装（通过有限差分计算雅可比）
    void assemble_matrix_general(Real dt,
                                const Vector& u,
                                const ComputationalDomain& domain) const;
    
    // 应用边界条件到矩阵和右端项
    void apply_boundary_to_system(Vector& rhs,
                                  const ComputationalDomain& domain) const;
};

} // namespace cfd
#endif // THETA_METHOD_H