#include "reconstructor/basic/SecondOrderReconstructor.h"
#include <algorithm>
#include <cmath>
#include <sstream>

namespace cfd {

void SecondOrderReconstructor::reconstruct(const Vector& q, 
                                         Vector& q_face_left,
                                         Vector& q_face_right,
                                         const ComputationalDomain& domain) {
    const int ist = domain.ist();
    const int ied = domain.ied();
    
    for (int i = ist; i <= ied; ++i) {
        int j = i - ist;
        
        // Central difference slope
        Real slope = 0.5 * (q[i] - q[i-2]);
        
        // Apply limiter parameter
        slope *= limiter_parameter_;
        
        // Apply slope limiter if enabled
        if (use_limiter_) {
            // Compute left and right slopes
            Real left_slope = q[i-1] - q[i-2];
            Real right_slope = q[i] - q[i-1];
            slope = minmod_limiter(left_slope, right_slope);
        }
        
        // Reconstruct face values
        q_face_left[j] = q[i-1] + 0.5 * slope;
        q_face_right[j] = q[i] - 0.5 * slope;
    }
}

std::string SecondOrderReconstructor::name() const { 
    std::ostringstream oss;
    oss << "Second-Order";
    
    if (use_limiter_) {
        oss << " (Minmod Limiter)";
    } else if (limiter_parameter_ == 1.0) {
        oss << " (Fromm)";
    } else if (limiter_parameter_ == 2.0) {
        oss << " (Beam-Warming)";
    } else if (limiter_parameter_ == 0.0) {
        oss << " (Upwind)";
    } else {
        oss << " (Param = " << limiter_parameter_ << ")";
    }
    
    return oss.str();
}

void SecondOrderReconstructor::apply_limiter(Vector& slopes) const {
    if (!use_limiter_) return;
    
    for (size_t i = 1; i < slopes.size() - 1; ++i) {
        Real left_slope = slopes[i] - slopes[i-1];
        Real right_slope = slopes[i+1] - slopes[i];
        slopes[i] = minmod_limiter(left_slope, right_slope);
    }
}

Real SecondOrderReconstructor::minmod_limiter(Real a, Real b) const {
    if (a * b <= 0.0) return 0.0;
    if (std::abs(a) < std::abs(b)) return a;
    return b;
}

Real SecondOrderReconstructor::van_leer_limiter(Real a, Real b) const {
    if (a * b <= 0.0) return 0.0;
    return 2.0 * a * b / (a + b);
}

Real SecondOrderReconstructor::superbee_limiter(Real a, Real b) const {
    if (a * b <= 0.0) return 0.0;
    Real r = a / b;
    Real phi = std::max(0.0, std::min(2.0 * r, 1.0));
    phi = std::max(phi, std::min(r, 2.0));
    return phi * b;
}

Real SecondOrderReconstructor::mc_limiter(Real a, Real b) const {
    if (a * b <= 0.0) return 0.0;
    Real c = (a + b) * 0.5;
    Real sign = (a > 0.0) ? 1.0 : -1.0;
    return sign * std::min(std::abs(c), std::min(2.0 * std::abs(a), 2.0 * std::abs(b)));
}

} // namespace cfd