#ifndef WENO_RECONSTRUCTOR_H
#define WENO_RECONSTRUCTOR_H

#include "../Reconstructor.h"

namespace cfd {

// ===================================================================
// WENO (Weighted Essentially Non-Oscillatory) reconstructor
// ===================================================================
class WenoReconstructor : public Reconstructor {
private:
    static constexpr Real eps_weno_ = 1.0e-6;
    int order_;
    int variant_;  // 0 = WENO-JS, 1 = WENO-Z, 2 = WENO-M
    
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
    
    // Get epsilon parameter
    static Real get_epsilon() { return eps_weno_; }
    
private:
    // WENO-JS weights
    Real weno_weight_js_L(Real v1, Real v2, Real v3) const;
    Real weno_weight_js_R(Real v1, Real v2, Real v3) const;
    
    // WENO-Z weights (improved version)
    Real weno_weight_z_L(Real v1, Real v2, Real v3) const;
    Real weno_weight_z_R(Real v1, Real v2, Real v3) const;
    
    // Smoothness indicators
    Real smoothness_indicator(Real v1, Real v2, Real v3) const;
    
    // Optimal weights (for WENO3)
    Real optimal_weight_L(int stencil) const;
    Real optimal_weight_R(int stencil) const;
};

} // namespace cfd

#endif // WENO_RECONSTRUCTOR_H