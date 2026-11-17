#ifndef ENO_RECONSTRUCTOR_H
#define ENO_RECONSTRUCTOR_H

#include "../Reconstructor.h"
#include <vector>

namespace cfd {

// ===================================================================
// ENO (Essentially Non-Oscillatory) reconstructor
// ===================================================================
class EnoReconstructor : public Reconstructor {
private:
    int spatial_order_;
    std::vector<int> lmc_;
    std::vector<std::vector<Real>> coef_;
    std::vector<std::vector<Real>> dd_;
    
    void initialize_coefficients(int order);
    
public:
    explicit EnoReconstructor(int order = 3);
    ~EnoReconstructor() override = default;
    
    void reconstruct(const Vector& q, 
                    Vector& q_face_left,
                    Vector& q_face_right,
                    const ComputationalDomain& domain) override;
    
    std::string name() const override;
    int order() const override { return spatial_order_; }
    int stencil_width() const override { return spatial_order_ + 1; }
    bool is_tvd() const override { return false; } // ENO is not strictly TVD
    bool is_monotonic() const override { return false; } // but essentially non-oscillatory
    bool is_linear() const override { return false; } // ENO is nonlinear due to stencil selection
};

} // namespace cfd

#endif // ENO_RECONSTRUCTOR_H