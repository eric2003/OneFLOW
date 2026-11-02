#ifndef COMPLEX_WAVE_INITIAL_CONDITION_H
#define COMPLEX_WAVE_INITIAL_CONDITION_H

#include "InitialCondition.h"
#include <cmath>
#include <string>

namespace cfd {

// ===================================================================
// 精确复杂波初始条件（源自文献的精确参数）
// 
// 参数说明：
//   a = 0.5       (区域4的中心)
//   z = -0.7      (区域1的中心)  
//   δ = 0.005     (高斯和椭圆的偏移量)
//   α = 10        (椭圆参数)
//   β = log(2) / (36δ²)  (高斯参数)
//
// 公式：
//   u_0(x) = (1/6)[G(x,β,z-δ) + G(x,β,z+δ) + 4G(x,β,z)],  -0.8 ≤ x ≤ -0.6
//            1,                                            -0.4 ≤ x ≤ -0.2
//            1 - |10(x-0.1)|,                               0 ≤ x ≤ 0.2
//            (1/6)[F(x,α,a-δ) + F(x,α,a+δ) + 4F(x,α,a)],   0.4 ≤ x ≤ 0.6
//            0,                                            其他区域
// 
// 其中：
//   G(x,β,z) = exp(-β(x-z)²)                     (高斯函数)
//   F(x,α,a) = √max(1 - α²(x-a)², 0)           (半椭圆函数)
// ===================================================================
class ComplexWaveInitialCondition : public InitialCondition {
private:
    // 精确参数
    struct PreciseParams {
        Real a = 0.5;          // 区域4的中心 (ellipse center)
        Real z = -0.7;         // 区域1的中心 (Gaussian center)
        Real delta = 0.005;    // 偏移量
        Real alpha = 10.0;     // 椭圆参数
        Real beta = 0.0;       // 高斯参数，将在构造函数中计算
        
        PreciseParams() {
            calculate_beta();
        }
        
        void calculate_beta() {
            static const Real LOG2 = std::log(2.0);  // 移除constexpr
            if (delta != 0.0) {
                beta = LOG2 / (36.0 * delta * delta);
            } else {
                beta = 0.0;
            }
        }
        
        void set_delta(Real new_delta) {
            delta = new_delta;
            calculate_beta();
        }
    };
    
    PreciseParams params_;
    std::string name_;
    bool initialized_;
    
    // 高斯函数: G(x,β,z) = exp(-β(x-z)²)
    Real gaussian(Real x, Real center) const {
        Real diff = x - center;
        return std::exp(-params_.beta * diff * diff);
    }
    
    // 半椭圆函数: F(x,α,a) = √max(1 - α²(x-a)², 0)
    Real half_ellipse(Real x, Real center) const {
        Real diff = x - center;
        Real arg = 1.0 - (params_.alpha * diff) * (params_.alpha * diff);
        return (arg > 0.0) ? std::sqrt(arg) : 0.0;
    }
    
public:
    // ===================================================================
    // 构造函数
    // ===================================================================
    
    // 默认构造函数（使用精确参数）
    ComplexWaveInitialCondition();
    
    // 自定义参数构造函数
    ComplexWaveInitialCondition(Real a, Real z, Real delta, 
                                Real alpha, const std::string& name = "Complex Wave");
    
    // ===================================================================
    // InitialCondition接口实现
    // ===================================================================
    
    void initialize(Vector& u_interior, const ComputationalDomain& domain) const override;
    
    Real evaluate(Real x) const override;
    
    std::string name() const override;
    
    int type_id() const override { return 4; }
    
    std::unique_ptr<InitialCondition> clone() const override;
    
    // ===================================================================
    // 参数访问器
    // ===================================================================
    
    Real a() const { return params_.a; }
    Real z() const { return params_.z; }
    Real delta() const { return params_.delta; }
    Real alpha() const { return params_.alpha; }
    Real beta() const { return params_.beta; }
    
    void set_a(Real val) { params_.a = val; }
    void set_z(Real val) { params_.z = val; }
    void set_delta(Real val) { params_.set_delta(val); }
    void set_alpha(Real val) { params_.alpha = val; }
    
    // ===================================================================
    // 辅助函数
    // ===================================================================
    
    // 打印参数信息
    void print_params() const;
    
    // 计算β值（静态方法）
    static Real calculate_beta(Real delta);
    
    // 创建不同分辨率的版本（用于收敛性研究）
    static std::unique_ptr<ComplexWaveInitialCondition> create_high_resolution();
    static std::unique_ptr<ComplexWaveInitialCondition> create_medium_resolution();
    static std::unique_ptr<ComplexWaveInitialCondition> create_low_resolution();
    
    // ===================================================================
    // 区域检查函数（用于调试）
    // ===================================================================
    
    // 判断点位于哪个区域
    enum Region {
        REGION_ZERO = 0,      // 其他区域
        REGION_GAUSSIAN = 1,  // 高斯组合区域
        REGION_SQUARE = 2,    // 方波区域
        REGION_TRIANGLE = 3,  // 三角波区域
        REGION_ELLIPSE = 4    // 椭圆组合区域
    };
    
    Region get_region(Real x) const;
    
    // 获取区域边界
    static constexpr Real region1_start() { return -0.8; }
    static constexpr Real region1_end()   { return -0.6; }
    static constexpr Real region2_start() { return -0.4; }
    static constexpr Real region2_end()   { return -0.2; }
    static constexpr Real region3_start() { return  0.0; }
    static constexpr Real region3_end()   { return  0.2; }
    static constexpr Real region4_start() { return  0.4; }
    static constexpr Real region4_end()   { return  0.6; }
    
private:
    // 初始化函数
    void initialize_params();
    
    // 区域计算函数
    Real calculate_region1(Real x) const;  // 高斯组合区域
    Real calculate_region2(Real x) const;  // 方波区域
    Real calculate_region3(Real x) const;  // 三角波区域
    Real calculate_region4(Real x) const;  // 椭圆组合区域
};

} // namespace cfd

#endif // COMPLEX_WAVE_INITIAL_CONDITION_H