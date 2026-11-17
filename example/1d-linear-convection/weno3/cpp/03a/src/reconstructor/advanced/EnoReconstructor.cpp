#include "reconstructor/advanced/EnoReconstructor.h"
#include <algorithm>
#include <stdexcept>
#include <sstream>

namespace cfd {

EnoReconstructor::EnoReconstructor(int order) : spatial_order_(order) {
    // 修改为支持1-7阶
    if (order < 1 || order > 7) {
        throw std::invalid_argument("ENO reconstructor supports order 1 through 7, got: " + std::to_string(order));
    }
    initialize_coefficients(order);
}

void EnoReconstructor::initialize_coefficients(int order) {
    // Initialize coefficient matrix
    // 注意：对于 order 阶重构，有 order+1 个模板，每个模板有 order 个系数
    coef_.resize(order + 1);
    for (auto& row : coef_) {
        row.resize(order, 0.0);
    }

    // Set coefficients based on order
    switch (order) {
    case 1:
        coef_[0][0] = 1.0;
        coef_[1][0] = 1.0;
        break;

    case 2:
        coef_[0][0] = 3.0/2.0;  coef_[0][1] = -1.0/2.0;
        coef_[1][0] = 1.0/2.0;  coef_[1][1] =  1.0/2.0;
        coef_[2][0] = -1.0/2.0; coef_[2][1] =  3.0/2.0;
        break;

    case 3:
        coef_[0][0] = 11.0/6.0; coef_[0][1] = -7.0/6.0; coef_[0][2] =  1.0/3.0;
        coef_[1][0] =  1.0/3.0; coef_[1][1] =  5.0/6.0; coef_[1][2] = -1.0/6.0;
        coef_[2][0] = -1.0/6.0; coef_[2][1] =  5.0/6.0; coef_[2][2] =  1.0/3.0;
        coef_[3][0] =  1.0/3.0; coef_[3][1] = -7.0/6.0; coef_[3][2] = 11.0/6.0;
        break;

    case 4:
        coef_[0][0] = 25.0/12.0; coef_[0][1] = -23.0/12.0; coef_[0][2] = 13.0/12.0; coef_[0][3] = -1.0/4.0;
        coef_[1][0] =  1.0/4.0;  coef_[1][1] =  13.0/12.0; coef_[1][2] = -5.0/12.0; coef_[1][3] =  1.0/12.0;
        coef_[2][0] = -1.0/12.0; coef_[2][1] =   7.0/12.0; coef_[2][2] =  7.0/12.0; coef_[2][3] = -1.0/12.0;
        coef_[3][0] =  1.0/12.0; coef_[3][1] =  -5.0/12.0; coef_[3][2] = 13.0/12.0; coef_[3][3] =  1.0/4.0;
        coef_[4][0] =  -1.0/4.0; coef_[4][1] =  13.0/12.0; coef_[4][2] = -23.0/12.0; coef_[4][3] = 25.0/12.0;
        break;

    case 5:
        coef_[0][0] = 137.0/60.0; coef_[0][1] = -163.0/60.0; coef_[0][2] = 137.0/60.0; coef_[0][3] = -21.0/20.0; coef_[0][4] =  1.0/5.0;
        coef_[1][0] =   1.0/5.0;  coef_[1][1] =   77.0/60.0; coef_[1][2] = -43.0/60.0; coef_[1][3] =  17.0/60.0; coef_[1][4] = -1.0/20.0;
        coef_[2][0] = -1.0/20.0;  coef_[2][1] =    9.0/20.0; coef_[2][2] =  47.0/60.0; coef_[2][3] = -13.0/60.0; coef_[2][4] =  1.0/30.0;
        coef_[3][0] =  1.0/30.0;  coef_[3][1] =  -13.0/60.0; coef_[3][2] =  47.0/60.0; coef_[3][3] =   9.0/20.0; coef_[3][4] = -1.0/20.0;
        coef_[4][0] = -1.0/20.0;  coef_[4][1] =   17.0/60.0; coef_[4][2] = -43.0/60.0; coef_[4][3] =  77.0/60.0; coef_[4][4] =   1.0/5.0;
        coef_[5][0] =   1.0/5.0;  coef_[5][1] =  -21.0/20.0; coef_[5][2] = 137.0/60.0; coef_[5][3] = -163.0/60.0; coef_[5][4] = 137.0/60.0;
        break;

    case 6:
        coef_[0][0] = 49.0/20.0;  coef_[0][1] = -71.0/20.0;  coef_[0][2] =  79.0/20.0;  coef_[0][3] = -163.0/60.0; coef_[0][4] =  31.0/30.0; coef_[0][5] = -1.0/6.0;
        coef_[1][0] =  1.0/6.0;   coef_[1][1] =  29.0/20.0;  coef_[1][2] = -21.0/20.0;  coef_[1][3] =  37.0/60.0;  coef_[1][4] = -13.0/60.0; coef_[1][5] =  1.0/30.0;
        coef_[2][0] = -1.0/30.0;  coef_[2][1] =  11.0/30.0;  coef_[2][2] =  19.0/20.0;  coef_[2][3] = -23.0/60.0;  coef_[2][4] =   7.0/60.0; coef_[2][5] = -1.0/60.0;
        coef_[3][0] =  1.0/60.0;  coef_[3][1] =  -2.0/15.0;  coef_[3][2] =  37.0/60.0;  coef_[3][3] =  37.0/60.0;  coef_[3][4] =  -2.0/15.0; coef_[3][5] =  1.0/60.0;
        coef_[4][0] = -1.0/60.0;  coef_[4][1] =   7.0/60.0;  coef_[4][2] = -23.0/60.0;  coef_[4][3] =  19.0/20.0;  coef_[4][4] =  11.0/30.0; coef_[4][5] = -1.0/30.0;
        coef_[5][0] =  1.0/30.0;  coef_[5][1] = -13.0/60.0;  coef_[5][2] =  37.0/60.0;  coef_[5][3] = -21.0/20.0;  coef_[5][4] =  29.0/20.0; coef_[5][5] =   1.0/6.0;
        coef_[6][0] = -1.0/6.0;   coef_[6][1] =  31.0/30.0;  coef_[6][2] = -163.0/60.0; coef_[6][3] =  79.0/20.0;  coef_[6][4] = -71.0/20.0; coef_[6][5] = 49.0/20.0;
        break;

    case 7:
        coef_[0][0] = 363.0/140.0; coef_[0][1] = -617.0/140.0; coef_[0][2] = 853.0/140.0; coef_[0][3] = -2341.0/420.0; coef_[0][4] = 667.0/210.0;  coef_[0][5] = -43.0/42.0;  coef_[0][6] = 1.0/7.0;
        coef_[1][0] =   1.0/7.0;   coef_[1][1] = 223.0/140.0;  coef_[1][2] = -197.0/140.0; coef_[1][3] = 153.0/140.0;   coef_[1][4] = -241.0/420.0; coef_[1][5] = 37.0/210.0;  coef_[1][6] = -1.0/42.0;
        coef_[2][0] =  -1.0/42.0;  coef_[2][1] =  13.0/42.0;   coef_[2][2] = 153.0/140.0;  coef_[2][3] = -241.0/420.0; coef_[2][4] = 109.0/420.0;  coef_[2][5] = -31.0/420.0; coef_[2][6] = 1.0/105.0;
        coef_[3][0] =  1.0/105.0;  coef_[3][1] = -19.0/210.0;  coef_[3][2] = 107.0/210.0;  coef_[3][3] = 319.0/420.0;   coef_[3][4] = -101.0/420.0; coef_[3][5] = 5.0/84.0;     coef_[3][6] = -1.0/140.0;
        coef_[4][0] = -1.0/140.0;  coef_[4][1] =   5.0/84.0;   coef_[4][2] = -101.0/420.0; coef_[4][3] = 319.0/420.0;   coef_[4][4] = 107.0/210.0;  coef_[4][5] = -19.0/210.0; coef_[4][6] = 1.0/105.0;
        coef_[5][0] =  1.0/105.0;  coef_[5][1] = -31.0/420.0;  coef_[5][2] = 109.0/420.0;  coef_[5][3] = -241.0/420.0; coef_[5][4] = 153.0/140.0;   coef_[5][5] = 13.0/42.0;    coef_[5][6] = -1.0/42.0;
        coef_[6][0] =  -1.0/42.0;  coef_[6][1] = 37.0/210.0;   coef_[6][2] = -241.0/420.0; coef_[6][3] = 153.0/140.0;   coef_[6][4] = -197.0/140.0; coef_[6][5] = 223.0/140.0;  coef_[6][6] = 1.0/7.0;
        coef_[7][0] =   1.0/7.0;   coef_[7][1] = -43.0/42.0;   coef_[7][2] = 667.0/210.0;  coef_[7][3] = -2341.0/420.0; coef_[7][4] = 853.0/140.0;  coef_[7][5] = -617.0/140.0; coef_[7][6] = 363.0/140.0;
        break;

    default:
        throw std::runtime_error("Unsupported ENO order: " + std::to_string(order));
    }
}

void EnoReconstructor::reconstruct(const Vector& q, 
                                 Vector& q_face_left,
                                 Vector& q_face_right,
                                 const ComputationalDomain& domain) {
    const int ntcells = domain.ntcells();
    const int ist = domain.ist();
    const int ied = domain.ied();
    
    // Resize working arrays if needed
    if (lmc_.size() != static_cast<size_t>(ntcells)) {
        lmc_.resize(ntcells, 0);
    }
    
    if (dd_.size() != static_cast<size_t>(spatial_order_)) {
        dd_.resize(spatial_order_);
        for (auto& row : dd_) {
            row.resize(ntcells, 0.0);
        }
    }
    
    // 1. Compute divided differences
    for (int i = 0; i < ntcells; ++i) {
        dd_[0][i] = q[i];
    }
    
    for (int m = 1; m < spatial_order_; ++m) {
        for (int j = 0; j < ntcells - m; ++j) {
            dd_[m][j] = dd_[m-1][j+1] - dd_[m-1][j];
        }
    }
    
    // 2. Select smoothest stencil
    for (int i = ist - 1; i <= ied; ++i) {
        lmc_[i] = i;
        for (int m = 1; m < spatial_order_; ++m) {
            if (std::abs(dd_[m][lmc_[i] - 1]) < std::abs(dd_[m][lmc_[i]])) {
                lmc_[i] = lmc_[i] - 1;
            }
        }
    }
    
    // 3. Reconstruct face values
    for (int i = ist; i <= ied; ++i) {
        int j = i - ist;
        
        int k1 = lmc_[i - 1];
        int k2 = lmc_[i];
        int r1 = (i - 1) - k1;
        int r2 = i - k2;
        
        q_face_left[j] = 0.0;
        q_face_right[j] = 0.0;
        
        for (int m = 0; m < spatial_order_; ++m) {
            q_face_left[j] += q[k1 + m] * coef_[r1 + 1][m];
            q_face_right[j] += q[k2 + m] * coef_[r2][m];
        }
    }
}

std::string EnoReconstructor::name() const { 
    std::ostringstream oss;
    oss << "ENO" << spatial_order_;
    return oss.str();
}

} // namespace cfd