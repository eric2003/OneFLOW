#include "initial_condition/ComplexWaveInitialCondition.h"
#include <iostream>
#include <sstream>
#include <iomanip>

namespace cfd {

// ===================================================================
// 构造函数
// ===================================================================

ComplexWaveInitialCondition::ComplexWaveInitialCondition() 
    : name_("Complex Wave (Precise)"), initialized_(false) {
    initialize_params();
}

ComplexWaveInitialCondition::ComplexWaveInitialCondition(Real a, Real z, Real delta, 
                                                         Real alpha, const std::string& name)
    : name_(name), initialized_(false) {
    params_.a = a;
    params_.z = z;
    params_.set_delta(delta);  // 这会自动计算beta
    params_.alpha = alpha;
    initialize_params();
}

void ComplexWaveInitialCondition::initialize_params() {
    if (!initialized_) {
        // 参数有效性检查
        if (params_.delta <= 0.0) {
            throw std::invalid_argument("Delta must be positive");
        }
        if (params_.alpha <= 0.0) {
            throw std::invalid_argument("Alpha must be positive");
        }
        
        initialized_ = true;
        print_params();  // 打印参数信息
    }
}

// ===================================================================
// 初始化函数
// ===================================================================

void ComplexWaveInitialCondition::initialize(Vector& u_interior, 
                                           const ComputationalDomain& domain) const {
    const auto& xcc = domain.mesh().cell_centers();
    
    for (int i = 0; i < domain.interior_cells(); ++i) {
        u_interior[i] = evaluate(xcc[i]);
    }
}

// ===================================================================
// 主评估函数
// ===================================================================

Real ComplexWaveInitialCondition::evaluate(Real x) const {
    // 区域1: 组合高斯脉冲 [-0.8, -0.6]
    if (x >= region1_start() && x <= region1_end()) {
        return calculate_region1(x);
    }
    
    // 区域2: 方波 [-0.4, -0.2]
    else if (x >= region2_start() && x <= region2_end()) {
        return calculate_region2(x);
    }
    
    // 区域3: 三角波 [0.0, 0.2]
    else if (x >= region3_start() && x <= region3_end()) {
        return calculate_region3(x);
    }
    
    // 区域4: 组合半椭圆 [0.4, 0.6]
    else if (x >= region4_start() && x <= region4_end()) {
        return calculate_region4(x);
    }
    
    // 其他区域: 零值
    else {
        return 0.0;
    }
}

// ===================================================================
// 区域计算函数
// ===================================================================

Real ComplexWaveInitialCondition::calculate_region1(Real x) const {
    // 区域1: 组合高斯脉冲
    // (1/6)[G(x,β,z-δ) + G(x,β,z+δ) + 4G(x,β,z)]
    Real g1 = gaussian(x, params_.z - params_.delta);  // G(x,β,z-δ)
    Real g2 = gaussian(x, params_.z + params_.delta);  // G(x,β,z+δ)
    Real g3 = gaussian(x, params_.z);                  // G(x,β,z)
    
    return (g1 + g2 + 4.0 * g3) / 6.0;
}

Real ComplexWaveInitialCondition::calculate_region2(Real x) const {
    // 区域2: 方波 (恒定值1)
    return 1.0;
}

Real ComplexWaveInitialCondition::calculate_region3(Real x) const {
    // 区域3: 三角波 1 - |10(x-0.1)|
    Real value = 1.0 - std::abs(10.0 * (x - 0.1));
    return (value > 0.0) ? value : 0.0;  // 处理浮点误差
}

Real ComplexWaveInitialCondition::calculate_region4(Real x) const {
    // 区域4: 组合半椭圆
    // (1/6)[F(x,α,a-δ) + F(x,α,a+δ) + 4F(x,α,a)]
    Real f1 = half_ellipse(x, params_.a - params_.delta);  // F(x,α,a-δ)
    Real f2 = half_ellipse(x, params_.a + params_.delta);  // F(x,α,a+δ)
    Real f3 = half_ellipse(x, params_.a);                  // F(x,α,a)
    
    return (f1 + f2 + 4.0 * f3) / 6.0;
}

// ===================================================================
// 其他成员函数实现
// ===================================================================

std::string ComplexWaveInitialCondition::name() const {
    return name_;
}

std::unique_ptr<InitialCondition> ComplexWaveInitialCondition::clone() const {
    return std::make_unique<ComplexWaveInitialCondition>(*this);
}

void ComplexWaveInitialCondition::print_params() const {
    std::cout << "=== Precise Complex Wave Parameters ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  a (ellipse center)     = " << params_.a << std::endl;
    std::cout << "  z (Gaussian center)    = " << params_.z << std::endl;
    std::cout << "  δ (offset)            = " << params_.delta << std::endl;
    std::cout << "  α (ellipse parameter) = " << params_.alpha << std::endl;
    std::cout << "  β (Gaussian parameter)= " << params_.beta 
              << " (log2/(36δ²))" << std::endl;
    std::cout << "  Region boundaries:" << std::endl;
    std::cout << "    Region 1 (Gaussian): [" << region1_start() << ", " << region1_end() << "]" << std::endl;
    std::cout << "    Region 2 (Square):   [" << region2_start() << ", " << region2_end() << "]" << std::endl;
    std::cout << "    Region 3 (Triangle): [" << region3_start() << ", " << region3_end() << "]" << std::endl;
    std::cout << "    Region 4 (Ellipse):  [" << region4_start() << ", " << region4_end() << "]" << std::endl;
    std::cout << std::endl;
}

Real ComplexWaveInitialCondition::calculate_beta(Real delta) {
    if (delta != 0.0) {
        return std::log(2.0) / (36.0 * delta * delta);
    }
    return 0.0;
}

ComplexWaveInitialCondition::Region ComplexWaveInitialCondition::get_region(Real x) const {
    if (x >= region1_start() && x <= region1_end()) {
        return REGION_GAUSSIAN;
    } else if (x >= region2_start() && x <= region2_end()) {
        return REGION_SQUARE;
    } else if (x >= region3_start() && x <= region3_end()) {
        return REGION_TRIANGLE;
    } else if (x >= region4_start() && x <= region4_end()) {
        return REGION_ELLIPSE;
    } else {
        return REGION_ZERO;
    }
}

// ===================================================================
// 工厂函数
// ===================================================================

std::unique_ptr<ComplexWaveInitialCondition> 
ComplexWaveInitialCondition::create_high_resolution() {
    // 高分辨率版本：更小的δ
    return std::make_unique<ComplexWaveInitialCondition>(
        0.5, -0.7, 0.001, 10.0, "Complex Wave (High Resolution)"
    );
}

std::unique_ptr<ComplexWaveInitialCondition> 
ComplexWaveInitialCondition::create_medium_resolution() {
    // 中等分辨率版本：默认参数
    return std::make_unique<ComplexWaveInitialCondition>();
}

std::unique_ptr<ComplexWaveInitialCondition> 
ComplexWaveInitialCondition::create_low_resolution() {
    // 低分辨率版本：更大的δ
    return std::make_unique<ComplexWaveInitialCondition>(
        0.5, -0.7, 0.01, 10.0, "Complex Wave (Low Resolution)"
    );
}

} // namespace cfd