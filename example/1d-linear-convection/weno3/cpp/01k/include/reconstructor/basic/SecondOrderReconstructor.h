#ifndef SECOND_ORDER_RECONSTRUCTOR_H
#define SECOND_ORDER_RECONSTRUCTOR_H

#include "../Reconstructor.h"

namespace cfd {

// ===================================================================
// Second-order linear reconstructor with optional limiter
// ===================================================================
class SecondOrderReconstructor : public Reconstructor {
private:
    Real limiter_parameter_;    // 1.0 = Fromm, 2.0 = Beam-Warming, 0.0 = Upwind
    bool use_limiter_;
    
public:
    SecondOrderReconstructor(Real limiter_param = 1.0, bool use_limiter = false)
        : limiter_parameter_(limiter_param), use_limiter_(use_limiter) {}
    
    ~SecondOrderReconstructor() override = default;
    
    void reconstruct(const Vector& q, 
                    Vector& q_face_left,
                    Vector& q_face_right,
                    const ComputationalDomain& domain) override;
    
    std::string name() const override;
    int order() const override { return 2; }
    int stencil_width() const override { return 3; }
    bool is_tvd() const override { return use_limiter_; }
    
    Real limiter_parameter() const { return limiter_parameter_; }
    bool use_limiter() const { return use_limiter_; }
    
    // Apply slope limiter (minmod, van Leer, etc.)
    void apply_limiter(Vector& slopes) const override;
    
private:
    // Different limiter functions
    Real minmod_limiter(Real a, Real b) const;
    Real van_leer_limiter(Real a, Real b) const;
    Real superbee_limiter(Real a, Real b) const;
    Real mc_limiter(Real a, Real b) const;  // monotonized central
};

} // namespace cfd

#endif // SECOND_ORDER_RECONSTRUCTOR_H