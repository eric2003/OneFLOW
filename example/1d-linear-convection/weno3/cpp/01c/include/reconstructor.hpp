// ==================== include/reconstructor.hpp ====================
#ifndef RECONSTRUCTOR_HPP
#define RECONSTRUCTOR_HPP

#include "domain.hpp"
#include "solution.hpp"
#include <vector>
#include <memory>
#include <string>
#include <cmath>

namespace cfd {

// Type aliases
using Real = double;
using Vector = std::vector<Real>;

// ===================================================================
// Reconstructor base class
// ===================================================================
class Reconstructor {
public:
    virtual ~Reconstructor() = default;
    
    // Reconstruct face values from cell averages
    virtual void reconstruct(const Vector& q, 
                           Vector& q_face_left,
                           Vector& q_face_right,
                           const ComputationalDomain& domain) = 0;
    
    // Get reconstructor name
    virtual std::string name() const = 0;
    
    // Get spatial order
    virtual int order() const = 0;
    
    // Check if reconstructor is linear
    virtual bool is_linear() const { return true; }
    
    // Check if reconstructor is nonlinear (e.g., WENO)
    virtual bool is_nonlinear() const { return !is_linear(); }
};

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
                    const ComputationalDomain& domain) override {
        const int ist = domain.ist();
        const int ied = domain.ied();
        
        for (int i = ist; i <= ied; ++i) {
            int j = i - ist;
            q_face_left[j] = q[i-1];
            q_face_right[j] = q[i];
        }
    }
    
    std::string name() const override { return "First-Order (Piecewise Constant)"; }
    int order() const override { return 1; }
};

// ===================================================================
// Second-order linear reconstructor
// ===================================================================
class SecondOrderReconstructor : public Reconstructor {
private:
    Real limiter_parameter_ = 1.0;  // 1.0 = Fromm, 2.0 = Beam-Warming, etc.
    
public:
    SecondOrderReconstructor() = default;
    explicit SecondOrderReconstructor(Real limiter_param) 
        : limiter_parameter_(limiter_param) {}
    
    void reconstruct(const Vector& q, 
                    Vector& q_face_left,
                    Vector& q_face_right,
                    const ComputationalDomain& domain) override {
        const int ist = domain.ist();
        const int ied = domain.ied();
        
        for (int i = ist; i <= ied; ++i) {
            int j = i - ist;
            
            // Central difference slope
            Real slope = 0.5 * (q[i] - q[i-2]);
            
            // Apply limiter parameter
            slope *= limiter_parameter_;
            
            // Reconstruct face values
            q_face_left[j] = q[i-1] + 0.5 * slope;
            q_face_right[j] = q[i] - 0.5 * slope;
        }
    }
    
    std::string name() const override { 
        if (limiter_parameter_ == 1.0) {
            return "Second-Order (Fromm)";
        } else {
            return "Second-Order (Param = " + std::to_string(limiter_parameter_) + ")";
        }
    }
    
    int order() const override { return 2; }
    Real limiter_parameter() const { return limiter_parameter_; }
};

// ===================================================================
// ENO reconstructor
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
    
    std::string name() const override { 
        return "ENO" + std::to_string(spatial_order_);
    }
    
    int order() const override { return spatial_order_; }
};

// ===================================================================
// WENO reconstructor
// ===================================================================
class WenoReconstructor : public Reconstructor {
private:
    static constexpr Real eps_weno_ = 1.0e-6;
    int order_;
    
    // WENO-JS weights for left-biased stencil
    Real weno_weight_L(Real v1, Real v2, Real v3) const;
    
    // WENO-JS weights for right-biased stencil
    Real weno_weight_R(Real v1, Real v2, Real v3) const;
    
public:
    explicit WenoReconstructor(int order = 3) : order_(order) {
        if (order != 3) {
            throw std::invalid_argument("Only WENO3 is currently implemented");
        }
    }
    
    ~WenoReconstructor() override = default;
    
    void reconstruct(const Vector& q, 
                    Vector& q_face_left,
                    Vector& q_face_right,
                    const ComputationalDomain& domain) override;
    
    std::string name() const override { 
        return "WENO" + std::to_string(order_);
    }
    
    int order() const override { return order_; }
    bool is_linear() const override { return false; }
    
    static Real get_epsilon() { return eps_weno_; }
};

// ===================================================================
// Reconstructor factory
// ===================================================================
class ReconstructorFactory {
public:
    static std::unique_ptr<Reconstructor> create_reconstructor(
        const std::string& scheme, int order = 3, Real param = 1.0) {
        
        if (scheme == "first" || scheme == "first-order") {
            return std::make_unique<FirstOrderReconstructor>();
        }
        else if (scheme == "second" || scheme == "second-order") {
            return std::make_unique<SecondOrderReconstructor>(param);
        }
        else if (scheme == "eno" || scheme == "ENO") {
            return std::make_unique<EnoReconstructor>(order);
        }
        else if (scheme == "weno" || scheme == "WENO") {
            return std::make_unique<WenoReconstructor>(order);
        }
        else {
            throw std::invalid_argument("Unknown reconstruction scheme: " + scheme);
        }
    }
    
    static std::unique_ptr<Reconstructor> create_reconstructor(
        const CfdConfig& config) {
        
        if (config.recon_scheme == "first" || config.recon_scheme == "first-order") {
            return std::make_unique<FirstOrderReconstructor>();
        }
        else if (config.recon_scheme == "second" || config.recon_scheme == "second-order") {
            return std::make_unique<SecondOrderReconstructor>(1.0);
        }
        else if (config.recon_scheme == "eno" || config.recon_scheme == "ENO") {
            return std::make_unique<EnoReconstructor>(config.spatial_order);
        }
        else if (config.recon_scheme == "weno" || config.recon_scheme == "WENO") {
            return std::make_unique<WenoReconstructor>(config.spatial_order);
        }
        else {
            throw std::invalid_argument("Unknown reconstruction scheme: " + config.recon_scheme);
        }
    }
    
    static std::vector<std::string> available_reconstructors() {
        return {
            "First-Order",
            "Second-Order",
            "ENO",
            "WENO"
        };
    }
};

} // namespace cfd

#endif // RECONSTRUCTOR_HPP