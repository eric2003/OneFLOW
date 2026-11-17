#ifndef FIRST_ORDER_RECONSTRUCTOR_H
#define FIRST_ORDER_RECONSTRUCTOR_H

#include "../Reconstructor.h"

namespace cfd {

// ===================================================================
// First-order reconstructor (piecewise constant)
// ===================================================================
class FirstOrderReconstructor : public Reconstructor {
public:
    FirstOrderReconstructor() = default;
    ~FirstOrderReconstructor() override = default;
    
    void reconstruct(const Vector& q, 
                    Vector& q_face_left,
                    Vector& q_face_right,
                    const ComputationalDomain& domain) override;
    
    std::string name() const override { return "First-Order (Piecewise Constant)"; }
    int order() const override { return 1; }
    int stencil_width() const override { return 1; }
    bool is_tvd() const override { return true; }
    bool is_monotonic() const override { return true; }
};

} // namespace cfd

#endif // FIRST_ORDER_RECONSTRUCTOR_H