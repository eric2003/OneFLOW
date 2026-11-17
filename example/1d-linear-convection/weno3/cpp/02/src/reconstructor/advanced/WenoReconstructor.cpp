#include "reconstructor/advanced/WenoReconstructor.h"
#include <cmath>
#include <sstream>

namespace cfd {

WenoReconstructor::WenoReconstructor(int order, int variant) 
    : order_(order), variant_(variant) {
    if (order != 3) {
        throw std::invalid_argument("Only WENO3 is currently implemented");
    }
}

void WenoReconstructor::reconstruct(const Vector& q, 
                                  Vector& q_face_left,
                                  Vector& q_face_right,
                                  const ComputationalDomain& domain) {
    const int ist = domain.ist();
    const int ied = domain.ied();
    
    // Left interface reconstruction
    for (int i = ist - 1; i < ied; ++i) {
        int j = i - (ist - 1);
        if (variant_ == 0) {
            q_face_left[j] = weno_weight_js_L(q[i-1], q[i], q[i+1]);
        } else {
            q_face_left[j] = weno_weight_z_L(q[i-1], q[i], q[i+1]);
        }
    }
    
    // Right interface reconstruction
    for (int i = ist; i <= ied; ++i) {
        int j = i - ist;
        if (variant_ == 0) {
            q_face_right[j] = weno_weight_js_R(q[i-1], q[i], q[i+1]);
        } else {
            q_face_right[j] = weno_weight_z_R(q[i-1], q[i], q[i+1]);
        }
    }
}

std::string WenoReconstructor::name() const { 
    std::ostringstream oss;
    oss << "WENO" << order_;
    
    switch (variant_) {
        case 0:
            oss << "-JS";
            break;
        case 1:
            oss << "-Z";
            break;
        case 2:
            oss << "-M";
            break;
    }
    
    return oss.str();
}

std::string WenoReconstructor::variant_name() const {
    switch (variant_) {
        case 0: return "WENO-JS (Jiang-Shu)";
        case 1: return "WENO-Z (Improved)";
        case 2: return "WENO-M (Mapped)";
        default: return "Unknown";
    }
}

Real WenoReconstructor::weno_weight_js_L(Real v1, Real v2, Real v3) const {
    // Smoothness indicators
    Real s0 = (v3 - v2) * (v3 - v2);
    Real s1 = (v2 - v1) * (v2 - v1);
    
    // Ideal weights (for WENO3)
    Real d0 = 2.0/3.0;
    Real d1 = 1.0/3.0;
    
    // Nonlinear weights
    Real alpha0 = d0 / ((eps_weno_ + s0) * (eps_weno_ + s0));
    Real alpha1 = d1 / ((eps_weno_ + s1) * (eps_weno_ + s1));
    
    Real sum_alpha = alpha0 + alpha1;
    Real w0 = alpha0 / sum_alpha;
    Real w1 = alpha1 / sum_alpha;
    
    // Candidate stencils
    Real q0 = 0.5 * v2 + 0.5 * v3;  // Right-biased
    Real q1 = -0.5 * v1 + 1.5 * v2; // Left-biased
    
    // Weighted combination
    return w0 * q0 + w1 * q1;
}

Real WenoReconstructor::weno_weight_js_R(Real v1, Real v2, Real v3) const {
    // Smoothness indicators
    Real s0 = (v2 - v1) * (v2 - v1);
    Real s1 = (v3 - v2) * (v3 - v2);
    
    // Ideal weights (for WENO3)
    Real d0 = 2.0/3.0;
    Real d1 = 1.0/3.0;
    
    // Nonlinear weights
    Real alpha0 = d0 / ((eps_weno_ + s0) * (eps_weno_ + s0));
    Real alpha1 = d1 / ((eps_weno_ + s1) * (eps_weno_ + s1));
    
    Real sum_alpha = alpha0 + alpha1;
    Real w0 = alpha0 / sum_alpha;
    Real w1 = alpha1 / sum_alpha;
    
    // Candidate stencils
    Real q0 = 0.5 * v1 + 0.5 * v2;  // Left-biased
    Real q1 = 1.5 * v2 - 0.5 * v3; // Right-biased
    
    // Weighted combination
    return w0 * q0 + w1 * q1;
}

Real WenoReconstructor::weno_weight_z_L(Real v1, Real v2, Real v3) const {
    // Smoothness indicators
    Real s0 = (v3 - v2) * (v3 - v2);
    Real s1 = (v2 - v1) * (v2 - v1);
    
    // Ideal weights
    Real d0 = 2.0/3.0;
    Real d1 = 1.0/3.0;
    
    // WENO-Z: use tau = |s0 - s1|
    Real tau = std::abs(s0 - s1);
    
    // Nonlinear weights
    Real alpha0 = d0 * (1.0 + (tau / (eps_weno_ + s0)));
    Real alpha1 = d1 * (1.0 + (tau / (eps_weno_ + s1)));
    
    Real sum_alpha = alpha0 + alpha1;
    Real w0 = alpha0 / sum_alpha;
    Real w1 = alpha1 / sum_alpha;
    
    // Candidate stencils
    Real q0 = 0.5 * v2 + 0.5 * v3;
    Real q1 = -0.5 * v1 + 1.5 * v2;
    
    return w0 * q0 + w1 * q1;
}

Real WenoReconstructor::weno_weight_z_R(Real v1, Real v2, Real v3) const {
    // Smoothness indicators
    Real s0 = (v2 - v1) * (v2 - v1);
    Real s1 = (v3 - v2) * (v3 - v2);
    
    // Ideal weights
    Real d0 = 2.0/3.0;
    Real d1 = 1.0/3.0;
    
    // WENO-Z: use tau = |s0 - s1|
    Real tau = std::abs(s0 - s1);
    
    // Nonlinear weights
    Real alpha0 = d0 * (1.0 + (tau / (eps_weno_ + s0)));
    Real alpha1 = d1 * (1.0 + (tau / (eps_weno_ + s1)));
    
    Real sum_alpha = alpha0 + alpha1;
    Real w0 = alpha0 / sum_alpha;
    Real w1 = alpha1 / sum_alpha;
    
    // Candidate stencils
    Real q0 = 0.5 * v1 + 0.5 * v2;
    Real q1 = 1.5 * v2 - 0.5 * v3;
    
    return w0 * q0 + w1 * q1;
}

Real WenoReconstructor::smoothness_indicator(Real v1, Real v2, Real v3) const {
    // For WENO3, the smoothness indicator is simple
    return (v3 - v2) * (v3 - v2) + (v2 - v1) * (v2 - v1);
}

Real WenoReconstructor::optimal_weight_L(int stencil) const {
    if (order_ == 3) {
        return (stencil == 0) ? 2.0/3.0 : 1.0/3.0;
    }
    return 1.0;
}

Real WenoReconstructor::optimal_weight_R(int stencil) const {
    if (order_ == 3) {
        return (stencil == 0) ? 1.0/3.0 : 2.0/3.0;
    }
    return 1.0;
}

} // namespace cfd