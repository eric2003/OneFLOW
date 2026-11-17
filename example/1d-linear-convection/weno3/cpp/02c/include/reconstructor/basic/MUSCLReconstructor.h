#ifndef MUSCL_RECONSTRUCTOR_H
#define MUSCL_RECONSTRUCTOR_H

#include "../Reconstructor.h"

namespace cfd {

    // ===================================================================
    // MUSCL (Monotonic Upstream-centered Scheme for Conservation Laws)
    // ===================================================================
    class MUSCLReconstructor : public Reconstructor {
    private:
        int order_;          // 2nd or 3rd order MUSCL
        std::string limiter_type_; // "minmod", "vanleer", "superbee", "mc"

    public:
        MUSCLReconstructor(int order = 2, const std::string& limiter = "minmod");
        ~MUSCLReconstructor() override = default;

        void reconstruct(const Vector& q, 
            Vector& q_face_left,
            Vector& q_face_right,
            const ComputationalDomain& domain) override;

        std::string name() const override;
        int order() const override { return order_; }
        int stencil_width() const override { return 3; }
        bool is_tvd() const override { return true; }
        bool is_monotonic() const override { return true; }
        bool is_linear() const override { return false; } // MUSCL with limiter is nonlinear

    private:
        // Limiter functions
        Real limiter(Real r) const;
        Real minmod(Real r) const;
        Real van_leer(Real r) const;
        Real superbee(Real r) const;
        Real mc(Real r) const;  // monotonized central

        // Compute limited slope (修正为3个参数)
        Real limited_slope(Real q_im1, Real q_i, Real q_ip1) const;
    };

} // namespace cfd

#endif // MUSCL_RECONSTRUCTOR_H