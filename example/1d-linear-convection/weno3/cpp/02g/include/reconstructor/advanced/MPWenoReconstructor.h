#ifndef MPWENO_RECONSTRUCTOR_H
#define MPWENO_RECONSTRUCTOR_H

#include "../advanced/WenoReconstructor.h"  // 修正包含路径

namespace cfd {

    // ===================================================================
    // MP-WENO (Monotonicity Preserving WENO)
    // ===================================================================
    class MPWenoReconstructor : public WenoReconstructor {
    private:
        Real mp_alpha_;  // MP parameter

    public:
        MPWenoReconstructor(int order = 3, Real mp_alpha = 2.0);
        ~MPWenoReconstructor() override = default;

        void reconstruct(const Vector& q, 
            Vector& q_face_left,
            Vector& q_face_right,
            const ComputationalDomain& domain) override;

        std::string name() const override;
        bool is_monotonic() const override { return true; }

    private:
        // MP-WENO specific methods
        Real mp_limiter(Real q_im2, Real q_im1, Real q_i, Real q_ip1, Real q_ip2) const;
        bool check_monotonicity(Real q_im1, Real q_i, Real q_ip1) const;
    };

} // namespace cfd

#endif // MPWENO_RECONSTRUCTOR_H