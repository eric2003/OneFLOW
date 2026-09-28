// ==================== src/reconstructor.cpp ====================
#include "reconstructor.hpp"
#include <algorithm>
#include <stdexcept>

namespace cfd {

// ===================================================================
// ENO reconstructor implementation
// ===================================================================
EnoReconstructor::EnoReconstructor(int order) : spatial_order_(order) {
    if (order < 1 || order > 3) {
        throw std::invalid_argument("ENO reconstructor only supports order 1, 2, or 3");
    }
    initialize_coefficients(order);
}

void EnoReconstructor::initialize_coefficients(int order) {
    // Initialize coefficient matrix
    coef_.resize(order + 1);
    for (auto& row : coef_) {
        row.resize(order, 0.0);
    }
    
    // Set coefficients based on order
    switch (order) {
        case 1:
            coef_[0][0] = 1.0;
            coef_[1][0] = 1.0;
            break;
            
        case 2:
            coef_[0][0] = 3.0/2.0;  coef_[0][1] = -1.0/2.0;
            coef_[1][0] = 1.0/2.0;  coef_[1][1] =  1.0/2.0;
            coef_[2][0] = -1.0/2.0; coef_[2][1] =  3.0/2.0;
            break;
            
        case 3:
            coef_[0][0] = 11.0/6.0; coef_[0][1] = -7.0/6.0; coef_[0][2] =  1.0/3.0;
            coef_[1][0] =  1.0/3.0; coef_[1][1] =  5.0/6.0; coef_[1][2] = -1.0/6.0;
            coef_[2][0] = -1.0/6.0; coef_[2][1] =  5.0/6.0; coef_[2][2] =  1.0/3.0;
            coef_[3][0] =  1.0/3.0; coef_[3][1] = -7.0/6.0; coef_[3][2] = 11.0/6.0;
            break;
    }
}

void EnoReconstructor::reconstruct(const Vector& q, 
                                  Vector& q_face_left,
                                  Vector& q_face_right,
                                  const ComputationalDomain& domain) {
    const int ntcells = domain.ntcells();
    const int ist = domain.ist();
    const int ied = domain.ied();
    
    // Resize working arrays if needed
    if (lmc_.size() != static_cast<size_t>(ntcells)) {
        lmc_.resize(ntcells, 0);
    }
    
    if (dd_.size() != static_cast<size_t>(spatial_order_)) {
        dd_.resize(spatial_order_);
        for (auto& row : dd_) {
            row.resize(ntcells, 0.0);
        }
    }
    
    // 1. Compute divided differences
    for (int i = 0; i < ntcells; ++i) {
        dd_[0][i] = q[i];
    }
    
    for (int m = 1; m < spatial_order_; ++m) {
        for (int j = 0; j < ntcells - m; ++j) {
            dd_[m][j] = dd_[m-1][j+1] - dd_[m-1][j];
        }
    }
    
    // 2. Select smoothest stencil
    for (int i = ist - 1; i <= ied; ++i) {
        lmc_[i] = i;
        for (int m = 1; m < spatial_order_; ++m) {
            if (std::abs(dd_[m][lmc_[i] - 1]) < std::abs(dd_[m][lmc_[i]])) {
                lmc_[i] = lmc_[i] - 1;
            }
        }
    }
    
    // 3. Reconstruct face values
    for (int i = ist; i <= ied; ++i) {
        int j = i - ist;
        
        int k1 = lmc_[i - 1];
        int k2 = lmc_[i];
        int r1 = (i - 1) - k1;
        int r2 = i - k2;
        
        q_face_left[j] = 0.0;
        q_face_right[j] = 0.0;
        
        for (int m = 0; m < spatial_order_; ++m) {
            q_face_left[j] += q[k1 + m] * coef_[r1 + 1][m];
            q_face_right[j] += q[k2 + m] * coef_[r2][m];
        }
    }
}

// ===================================================================
// WENO reconstructor implementation
// ===================================================================
Real WenoReconstructor::weno_weight_L(Real v1, Real v2, Real v3) const {
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

Real WenoReconstructor::weno_weight_R(Real v1, Real v2, Real v3) const {
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

void WenoReconstructor::reconstruct(const Vector& q, 
                                   Vector& q_face_left,
                                   Vector& q_face_right,
                                   const ComputationalDomain& domain) {
    const int ist = domain.ist();
    const int ied = domain.ied();
    
    // Left interface reconstruction
    for (int i = ist - 1; i < ied; ++i) {
        int j = i - (ist - 1);
        q_face_left[j] = weno_weight_L(q[i-1], q[i], q[i+1]);
    }
    
    // Right interface reconstruction
    for (int i = ist; i <= ied; ++i) {
        int j = i - ist;
        q_face_right[j] = weno_weight_R(q[i-1], q[i], q[i+1]);
    }
}

} // namespace cfd