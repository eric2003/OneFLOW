#ifndef WENO_RECONSTRUCTOR_H
#define WENO_RECONSTRUCTOR_H

#include "../Reconstructor.h"
#include <array>
#include <vector>

namespace cfd {

    // ===================================================================
    // WENO (Weighted Essentially Non-Oscillatory) reconstructor
    // 支持 3阶 和 5阶 WENO
    // ===================================================================
    class WenoReconstructor : public Reconstructor {
    private:
        static constexpr Real eps_weno_ = 1.0e-40;  // WENO5需要更小的epsilon
        int order_;          // 3 或 5
        int variant_;        // 0 = WENO-JS, 1 = WENO-Z

    public:
        explicit WenoReconstructor(int order = 3, int variant = 0);
        ~WenoReconstructor() override = default;

        void reconstruct(const Vector& q, 
            Vector& q_face_left,
            Vector& q_face_right,
            const ComputationalDomain& domain) override;

        std::string name() const override;
        int order() const override { return order_; }
        int stencil_width() const override { return order_ + 1; }
        bool is_tvd() const override { return false; }
        bool is_monotonic() const override { return false; }
        bool is_linear() const override { return false; }

        // WENO variant name
        std::string variant_name() const;

        // Check if order is supported
        static bool is_order_supported(int order) {
            return (order == 3 || order == 5);
        }

        // Get epsilon parameter
        static Real get_epsilon() { return eps_weno_; }

    private:
        // WENO3 reconstruction
        void reconstruct_weno3(const Vector& q, 
            Vector& q_face_left,
            Vector& q_face_right,
            const ComputationalDomain& domain) const;

        // WENO5 reconstruction
        void reconstruct_weno5(const Vector& q, 
            Vector& q_face_left,
            Vector& q_face_right,
            const ComputationalDomain& domain) const;

        // ===================================================================
        // WENO3 specific methods
        // ===================================================================

        // WENO3-JS weights
        Real weno3_weight_js_L(Real v1, Real v2, Real v3) const;
        Real weno3_weight_js_R(Real v1, Real v2, Real v3) const;

        // WENO3-Z weights
        Real weno3_weight_z_L(Real v1, Real v2, Real v3) const;
        Real weno3_weight_z_R(Real v1, Real v2, Real v3) const;

        // WENO3 smoothness indicator
        Real smoothness_indicator_3(Real v1, Real v2, Real v3) const;

        // WENO3 optimal weights
        Real optimal_weight_3_L(int stencil) const;
        Real optimal_weight_3_R(int stencil) const;

        // ===================================================================
        // WENO5 specific methods
        // ===================================================================

        // WENO5-JS weights
        std::array<Real, 3> weno5_weights_js_L(Real v1, Real v2, Real v3, Real v4, Real v5) const;
        std::array<Real, 3> weno5_weights_js_R(Real v1, Real v2, Real v3, Real v4, Real v5) const;

        // WENO5-Z weights
        std::array<Real, 3> weno5_weights_z_L(Real v1, Real v2, Real v3, Real v4, Real v5) const;
        std::array<Real, 3> weno5_weights_z_R(Real v1, Real v2, Real v3, Real v4, Real v5) const;

        // WENO5 smoothness indicators (3 stencils)
        std::array<Real, 3> smoothness_indicators_5(Real v1, Real v2, Real v3, Real v4, Real v5) const;

        // WENO5 optimal weights
        static constexpr std::array<Real, 3> optimal_weights_5_L = {0.1, 0.6, 0.3};
        static constexpr std::array<Real, 3> optimal_weights_5_R = {0.3, 0.6, 0.1};

        // WENO5 candidate stencils for left face
        std::array<Real, 3> candidate_stencils_5_L(Real v1, Real v2, Real v3, Real v4, Real v5) const;

        // WENO5 candidate stencils for right face
        std::array<Real, 3> candidate_stencils_5_R(Real v1, Real v2, Real v3, Real v4, Real v5) const;

        // ===================================================================
        // Common helper methods
        // ===================================================================

        // Get epsilon with safeguard
        Real get_epsilon_safe() const {
            return std::max(eps_weno_, std::numeric_limits<Real>::epsilon() * 10.0);
        }

        // Power function for weights
        Real power_for_weights() const {
            return (order_ == 3) ? 2.0 : 4.0;  // WENO3用2次方，WENO5用4次方
        }
    };

} // namespace cfd

#endif // WENO_RECONSTRUCTOR_H