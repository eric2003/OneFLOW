#include "reconstructor/advanced/WenoReconstructor.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace cfd {

    // ===================================================================
    // Constructor
    // ===================================================================

    WenoReconstructor::WenoReconstructor(int order, int variant) 
        : order_(order), variant_(variant) {
        if (!is_order_supported(order_)) {
            throw std::invalid_argument("WENO reconstructor only supports order 3 or 5");
        }
        if (variant_ < 0 || variant_ > 1) {
            throw std::invalid_argument("WENO variant must be 0 (JS) or 1 (Z)");
        }
    }

    // ===================================================================
    // Main reconstruction method
    // ===================================================================

    void WenoReconstructor::reconstruct(const Vector& q, 
        Vector& q_face_left,
        Vector& q_face_right,
        const ComputationalDomain& domain) {
        if (order_ == 3) {
            reconstruct_weno3(q, q_face_left, q_face_right, domain);
        } else if (order_ == 5) {
            reconstruct_weno5(q, q_face_left, q_face_right, domain);
        } else {
            throw std::runtime_error("Unsupported WENO order: " + std::to_string(order_));
        }
    }

    // ===================================================================
    // WENO3 reconstruction
    // ===================================================================

    void WenoReconstructor::reconstruct_weno3(const Vector& q, 
        Vector& q_face_left,
        Vector& q_face_right,
        const ComputationalDomain& domain) const {
        const int ist = domain.ist();
        const int ied = domain.ied();

        // Left interface reconstruction
        for (int i = ist - 1; i < ied; ++i) {
            int j = i - (ist - 1);
            if (variant_ == 0) {
                q_face_left[j] = weno3_weight_js_L(q[i-1], q[i], q[i+1]);
            } else {
                q_face_left[j] = weno3_weight_z_L(q[i-1], q[i], q[i+1]);
            }
        }

        // Right interface reconstruction
        for (int i = ist; i <= ied; ++i) {
            int j = i - ist;
            if (variant_ == 0) {
                q_face_right[j] = weno3_weight_js_R(q[i-1], q[i], q[i+1]);
            } else {
                q_face_right[j] = weno3_weight_z_R(q[i-1], q[i], q[i+1]);
            }
        }
    }

    // ===================================================================
    // WENO5 reconstruction
    // ===================================================================

    void WenoReconstructor::reconstruct_weno5(const Vector& q, 
        Vector& q_face_left,
        Vector& q_face_right,
        const ComputationalDomain& domain) const {
        const int ist = domain.ist();
        const int ied = domain.ied();

        // Left interface reconstruction
        for (int i = ist - 1; i < ied; ++i) {
            int j = i - (ist - 1);

            // Get values for 5-point stencil: i-2, i-1, i, i+1, i+2
            Real v1 = q[i-2];
            Real v2 = q[i-1];
            Real v3 = q[i];
            Real v4 = q[i+1];
            Real v5 = q[i+2];

            std::array<Real, 3> weights;
            if (variant_ == 0) {
                weights = weno5_weights_js_L(v1, v2, v3, v4, v5);
            } else {
                weights = weno5_weights_z_L(v1, v2, v3, v4, v5);
            }

            auto candidates = candidate_stencils_5_L(v1, v2, v3, v4, v5);
            q_face_left[j] = weights[0] * candidates[0] + 
                weights[1] * candidates[1] + 
                weights[2] * candidates[2];
        }

        // Right interface reconstruction
        for (int i = ist; i <= ied; ++i) {
            int j = i - ist;

            // Get values for 5-point stencil: i-2, i-1, i, i+1, i+2
            Real v1 = q[i-2];
            Real v2 = q[i-1];
            Real v3 = q[i];
            Real v4 = q[i+1];
            Real v5 = q[i+2];

            std::array<Real, 3> weights;
            if (variant_ == 0) {
                weights = weno5_weights_js_R(v1, v2, v3, v4, v5);
            } else {
                weights = weno5_weights_z_R(v1, v2, v3, v4, v5);
            }

            auto candidates = candidate_stencils_5_R(v1, v2, v3, v4, v5);
            q_face_right[j] = weights[0] * candidates[0] + 
                weights[1] * candidates[1] + 
                weights[2] * candidates[2];
        }
    }

    // ===================================================================
    // WENO3 specific methods
    // ===================================================================

    Real WenoReconstructor::weno3_weight_js_L(Real v1, Real v2, Real v3) const {
        // Smoothness indicators
        Real s0 = (v3 - v2) * (v3 - v2);
        Real s1 = (v2 - v1) * (v2 - v1);

        // Ideal weights
        Real d0 = 2.0/3.0;
        Real d1 = 1.0/3.0;

        // Nonlinear weights
        Real eps = get_epsilon_safe();
        Real alpha0 = d0 / ((eps + s0) * (eps + s0));
        Real alpha1 = d1 / ((eps + s1) * (eps + s1));

        Real sum_alpha = alpha0 + alpha1;
        Real w0 = alpha0 / sum_alpha;
        Real w1 = alpha1 / sum_alpha;

        // Candidate stencils
        Real q0 = 0.5 * v2 + 0.5 * v3;  // Right-biased
        Real q1 = -0.5 * v1 + 1.5 * v2; // Left-biased

        return w0 * q0 + w1 * q1;
    }

    Real WenoReconstructor::weno3_weight_js_R(Real v1, Real v2, Real v3) const {
        // Smoothness indicators
        Real s0 = (v2 - v1) * (v2 - v1);
        Real s1 = (v3 - v2) * (v3 - v2);

        // Ideal weights
        Real d0 = 2.0/3.0;
        Real d1 = 1.0/3.0;

        // Nonlinear weights
        Real eps = get_epsilon_safe();
        Real alpha0 = d0 / ((eps + s0) * (eps + s0));
        Real alpha1 = d1 / ((eps + s1) * (eps + s1));

        Real sum_alpha = alpha0 + alpha1;
        Real w0 = alpha0 / sum_alpha;
        Real w1 = alpha1 / sum_alpha;

        // Candidate stencils
        Real q0 = 0.5 * v1 + 0.5 * v2;  // Left-biased
        Real q1 = 1.5 * v2 - 0.5 * v3;  // Right-biased

        return w0 * q0 + w1 * q1;
    }

    Real WenoReconstructor::weno3_weight_z_L(Real v1, Real v2, Real v3) const {
        // Smoothness indicators
        Real s0 = (v3 - v2) * (v3 - v2);
        Real s1 = (v2 - v1) * (v2 - v1);

        // Ideal weights
        Real d0 = 2.0/3.0;
        Real d1 = 1.0/3.0;

        // WENO-Z parameter
        Real tau = std::abs(s0 - s1);

        // Nonlinear weights
        Real eps = get_epsilon_safe();
        Real alpha0 = d0 * (1.0 + (tau / (eps + s0)));
        Real alpha1 = d1 * (1.0 + (tau / (eps + s1)));

        Real sum_alpha = alpha0 + alpha1;
        Real w0 = alpha0 / sum_alpha;
        Real w1 = alpha1 / sum_alpha;

        // Candidate stencils
        Real q0 = 0.5 * v2 + 0.5 * v3;
        Real q1 = -0.5 * v1 + 1.5 * v2;

        return w0 * q0 + w1 * q1;
    }

    Real WenoReconstructor::weno3_weight_z_R(Real v1, Real v2, Real v3) const {
        // Smoothness indicators
        Real s0 = (v2 - v1) * (v2 - v1);
        Real s1 = (v3 - v2) * (v3 - v2);

        // Ideal weights
        Real d0 = 2.0/3.0;
        Real d1 = 1.0/3.0;

        // WENO-Z parameter
        Real tau = std::abs(s0 - s1);

        // Nonlinear weights
        Real eps = get_epsilon_safe();
        Real alpha0 = d0 * (1.0 + (tau / (eps + s0)));
        Real alpha1 = d1 * (1.0 + (tau / (eps + s1)));

        Real sum_alpha = alpha0 + alpha1;
        Real w0 = alpha0 / sum_alpha;
        Real w1 = alpha1 / sum_alpha;

        // Candidate stencils
        Real q0 = 0.5 * v1 + 0.5 * v2;
        Real q1 = 1.5 * v2 - 0.5 * v3;

        return w0 * q0 + w1 * q1;
    }

    Real WenoReconstructor::smoothness_indicator_3(Real v1, Real v2, Real v3) const {
        return (v3 - v2) * (v3 - v2) + (v2 - v1) * (v2 - v1);
    }

    Real WenoReconstructor::optimal_weight_3_L(int stencil) const {
        return (stencil == 0) ? 2.0/3.0 : 1.0/3.0;
    }

    Real WenoReconstructor::optimal_weight_3_R(int stencil) const {
        return (stencil == 0) ? 1.0/3.0 : 2.0/3.0;
    }

    // ===================================================================
    // WENO5 specific methods
    // ===================================================================

    std::array<Real, 3> WenoReconstructor::weno5_weights_js_L(Real v1, Real v2, Real v3, 
        Real v4, Real v5) const {
        auto beta = smoothness_indicators_5(v1, v2, v3, v4, v5);
        Real eps = get_epsilon_safe();
        Real p = power_for_weights();

        std::array<Real, 3> alpha;
        for (int k = 0; k < 3; ++k) {
            alpha[k] = optimal_weights_5_L[k] / std::pow(eps + beta[k], p);
        }

        Real sum_alpha = alpha[0] + alpha[1] + alpha[2];
        std::array<Real, 3> weights;
        for (int k = 0; k < 3; ++k) {
            weights[k] = alpha[k] / sum_alpha;
        }

        return weights;
    }

    std::array<Real, 3> WenoReconstructor::weno5_weights_js_R(Real v1, Real v2, Real v3, 
        Real v4, Real v5) const {
        auto beta = smoothness_indicators_5(v1, v2, v3, v4, v5);
        Real eps = get_epsilon_safe();
        Real p = power_for_weights();

        std::array<Real, 3> alpha;
        for (int k = 0; k < 3; ++k) {
            alpha[k] = optimal_weights_5_R[k] / std::pow(eps + beta[k], p);
        }

        Real sum_alpha = alpha[0] + alpha[1] + alpha[2];
        std::array<Real, 3> weights;
        for (int k = 0; k < 3; ++k) {
            weights[k] = alpha[k] / sum_alpha;
        }

        return weights;
    }

    std::array<Real, 3> WenoReconstructor::weno5_weights_z_L(Real v1, Real v2, Real v3, 
        Real v4, Real v5) const {
        auto beta = smoothness_indicators_5(v1, v2, v3, v4, v5);
        Real eps = get_epsilon_safe();

        // WENO-Z parameter: tau5 = |¦Â0 - ¦Â2|
        Real tau5 = std::abs(beta[0] - beta[2]);

        std::array<Real, 3> alpha;
        for (int k = 0; k < 3; ++k) {
            alpha[k] = optimal_weights_5_L[k] * (1.0 + (tau5 / (eps + beta[k])));
        }

        Real sum_alpha = alpha[0] + alpha[1] + alpha[2];
        std::array<Real, 3> weights;
        for (int k = 0; k < 3; ++k) {
            weights[k] = alpha[k] / sum_alpha;
        }

        return weights;
    }

    std::array<Real, 3> WenoReconstructor::weno5_weights_z_R(Real v1, Real v2, Real v3, 
        Real v4, Real v5) const {
        auto beta = smoothness_indicators_5(v1, v2, v3, v4, v5);
        Real eps = get_epsilon_safe();

        // WENO-Z parameter: tau5 = |¦Â0 - ¦Â2|
        Real tau5 = std::abs(beta[0] - beta[2]);

        std::array<Real, 3> alpha;
        for (int k = 0; k < 3; ++k) {
            alpha[k] = optimal_weights_5_R[k] * (1.0 + (tau5 / (eps + beta[k])));
        }

        Real sum_alpha = alpha[0] + alpha[1] + alpha[2];
        std::array<Real, 3> weights;
        for (int k = 0; k < 3; ++k) {
            weights[k] = alpha[k] / sum_alpha;
        }

        return weights;
    }

    std::array<Real, 3> WenoReconstructor::smoothness_indicators_5(Real v1, Real v2, Real v3, 
        Real v4, Real v5) const {
        std::array<Real, 3> beta;

        // ¦Â0 for stencil S0 = {i-2, i-1, i}
        beta[0] = (13.0/12.0) * std::pow(v1 - 2.0*v2 + v3, 2) 
            + (1.0/4.0) * std::pow(v1 - 4.0*v2 + 3.0*v3, 2);

        // ¦Â1 for stencil S1 = {i-1, i, i+1}
        beta[1] = (13.0/12.0) * std::pow(v2 - 2.0*v3 + v4, 2) 
            + (1.0/4.0) * std::pow(v2 - v4, 2);

        // ¦Â2 for stencil S2 = {i, i+1, i+2}
        beta[2] = (13.0/12.0) * std::pow(v3 - 2.0*v4 + v5, 2) 
            + (1.0/4.0) * std::pow(3.0*v3 - 4.0*v4 + v5, 2);

        return beta;
    }

    std::array<Real, 3> WenoReconstructor::candidate_stencils_5_L(Real v1, Real v2, Real v3, 
        Real v4, Real v5) const {
        std::array<Real, 3> candidates;

        // Candidate stencils for left face reconstruction
        // S0: {i-2, i-1, i} -> q0 = (1/3)v1 - (7/6)v2 + (11/6)v3
        candidates[0] = (1.0/3.0)*v1 - (7.0/6.0)*v2 + (11.0/6.0)*v3;

        // S1: {i-1, i, i+1} -> q1 = -(1/6)v2 + (5/6)v3 + (1/3)v4
        candidates[1] = -(1.0/6.0)*v2 + (5.0/6.0)*v3 + (1.0/3.0)*v4;

        // S2: {i, i+1, i+2} -> q2 = (1/3)v3 + (5/6)v4 - (1/6)v5
        candidates[2] = (1.0/3.0)*v3 + (5.0/6.0)*v4 - (1.0/6.0)*v5;

        return candidates;
    }

    std::array<Real, 3> WenoReconstructor::candidate_stencils_5_R(Real v1, Real v2, Real v3, 
        Real v4, Real v5) const {
        std::array<Real, 3> candidates;

        // Candidate stencils for right face reconstruction
        // S0: {i-2, i-1, i} -> q0 = (1/3)v1 + (5/6)v2 - (1/6)v3
        candidates[0] = (1.0/3.0)*v1 + (5.0/6.0)*v2 - (1.0/6.0)*v3;

        // S1: {i-1, i, i+1} -> q1 = -(1/6)v2 + (5/6)v3 + (1/3)v4
        candidates[1] = -(1.0/6.0)*v2 + (5.0/6.0)*v3 + (1.0/3.0)*v4;

        // S2: {i, i+1, i+2} -> q2 = (1/3)v3 - (7/6)v4 + (11/6)v5
        candidates[2] = (1.0/3.0)*v3 - (7.0/6.0)*v4 + (11.0/6.0)*v5;

        return candidates;
    }

    // ===================================================================
    // Name and variant name
    // ===================================================================

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
        default:
            oss << " (unknown variant)";
        }

        return oss.str();
    }

    std::string WenoReconstructor::variant_name() const {
        switch (variant_) {
        case 0: return "WENO-JS (Jiang-Shu)";
        case 1: return "WENO-Z (Improved)";
        default: return "Unknown";
        }
    }

} // namespace cfd