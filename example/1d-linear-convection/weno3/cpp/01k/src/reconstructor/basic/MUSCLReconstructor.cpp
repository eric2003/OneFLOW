#include "reconstructor/basic/MUSCLReconstructor.h"
#include <cmath>
#include <algorithm>
#include <sstream>

namespace cfd {

    MUSCLReconstructor::MUSCLReconstructor(int order, const std::string& limiter)
        : order_(order), limiter_type_(limiter) {
        if (order != 2 && order != 3) {
            throw std::invalid_argument("MUSCL reconstructor only supports order 2 or 3");
        }
    }

    void MUSCLReconstructor::reconstruct(const Vector& q, 
        Vector& q_face_left,
        Vector& q_face_right,
        const ComputationalDomain& domain) {
        const int ist = domain.ist();
        const int ied = domain.ied();

        for (int i = ist; i <= ied; ++i) {
            int j = i - ist;

            // Get neighboring values
            Real q_im2 = q[i-2];
            Real q_im1 = q[i-1];
            Real q_i   = q[i];
            Real q_ip1 = (i < ied) ? q[i+1] : q[i];

            if (order_ == 2) {
                // 2nd order MUSCL - use slope limiter
                Real slope_left = q_i - q_im1;
                Real slope_right = q_im1 - q_im2;

                // Apply limiter using slope ratio
                Real r_left = (std::abs(slope_right) > 1e-12) ? slope_left / slope_right : 1.0;
                Real phi_left = limiter(r_left);
                Real limited_slope_left = phi_left * slope_left;

                // Reconstruct face values
                q_face_left[j] = q_im1 + 0.5 * limited_slope_left;
                q_face_right[j] = q_i - 0.5 * limited_slope_left;
            } 
            else {
                // 3rd order MUSCL - use limited_slope function
                Real slope = limited_slope(q_im1, q_i, q_ip1);

                q_face_left[j] = q_im1 + 0.5 * slope;
                q_face_right[j] = q_i - 0.5 * slope;
            }
        }
    }

    std::string MUSCLReconstructor::name() const { 
        std::ostringstream oss;
        oss << "MUSCL" << order_ << " (" << limiter_type_ << " limiter)";
        return oss.str();
    }

    Real MUSCLReconstructor::limiter(Real r) const {
        if (limiter_type_ == "minmod") {
            return minmod(r);
        } else if (limiter_type_ == "vanleer") {
            return van_leer(r);
        } else if (limiter_type_ == "superbee") {
            return superbee(r);
        } else if (limiter_type_ == "mc") {
            return mc(r);
        } else {
            // Default to minmod
            return minmod(r);
        }
    }

    Real MUSCLReconstructor::minmod(Real r) const {
        return std::max(0.0, std::min(1.0, r));
    }

    Real MUSCLReconstructor::van_leer(Real r) const {
        return (r + std::abs(r)) / (1.0 + std::abs(r));
    }

    Real MUSCLReconstructor::superbee(Real r) const {
        // superbee returns max(0, min(2r, 1), min(r, 2))
        Real term1 = std::max(0.0, std::min(2.0 * r, 1.0));
        return std::max(term1, std::min(r, 2.0));
    }

    //Real MUSCLReconstructor::mc(Real r) const {
    //    return std::max(0.0, std::min(2.0 * r, 0.5 * (1.0 + r), 2.0));
    //}

    Real MUSCLReconstructor::mc(Real r) const {
        return std::max(0.0, std::min({2.0 * r, 0.5 * (1.0 + r), 2.0}));
    }

    Real MUSCLReconstructor::limited_slope(Real q_im1, Real q_i, Real q_ip1) const {
        // Compute three candidate slopes for MUSCL3
        Real slope_left = q_i - q_im1;
        Real slope_center = 0.5 * (q_ip1 - q_im1);
        Real slope_right = q_ip1 - q_i;

        // Apply minmod limiter to ensure TVD property
        if (slope_left * slope_center <= 0.0 || slope_left * slope_right <= 0.0) {
            return 0.0;
        }

        Real sign = (slope_left > 0.0) ? 1.0 : -1.0;

        if (limiter_type_ == "minmod") {
            // Minmod: take the smallest absolute value with the same sign
            Real abs1 = std::abs(slope_left);
            Real abs2 = std::abs(slope_center);
            Real abs3 = std::abs(slope_right);
            return sign * std::min(std::min(abs1, abs2), abs3);
        } 
        else if (limiter_type_ == "vanleer") {
            // Van Leer harmonic mean limiter
            Real s12 = slope_left * slope_center;
            Real s13 = slope_left * slope_right;
            return sign * (2.0 * s12 * s13) / (s12 + s13 + 1e-12);
        }
        else if (limiter_type_ == "superbee") {
            // Superbee limiter
            Real abs1 = std::abs(slope_left);
            Real abs2 = std::abs(slope_center);
            Real abs3 = std::abs(slope_right);
            Real min_abs = std::min(std::min(abs1, abs2), abs3);
            Real max_abs = std::max(std::max(abs1, abs2), abs3);
            return sign * (min_abs + max_abs) * 0.5;
        }
        else { 
            // mc (monotonized central) or default
            Real abs1 = std::abs(slope_left);
            Real abs2 = std::abs(slope_center);
            Real abs3 = std::abs(slope_right);
            return sign * std::min(abs1, std::min(2.0 * abs2, 2.0 * abs3));
        }
    }

} // namespace cfd