#ifndef INITIALIZER_HPP
#define INITIALIZER_HPP

#include "domain.hpp"
#include <vector>
#include <functional>
#include <cmath>
#include <memory>
#include <string>

namespace cfd {

// Type aliases
using Real = double;
using Vector = std::vector<Real>;

// ===================================================================
// Initial condition base class
// ===================================================================
class InitialCondition {
public:
    virtual ~InitialCondition() = default;
    
    // Initialize solution field at cell centers
    virtual void initialize(Vector& u_interior, const ComputationalDomain& domain) const = 0;
    
    // Initialize solution field with ghost cells
    virtual void initialize_with_ghosts(Vector& u_with_ghosts, const ComputationalDomain& domain) const {
        // Get interior solution
        Vector u_interior = domain.extract_interior(u_with_ghosts);
        
        // Initialize interior
        initialize(u_interior, domain);
        
        // Copy back to ghosted array
        const int ist = domain.ist();
        for (int i = 0; i < domain.interior_cells(); ++i) {
            u_with_ghosts[ist + i] = u_interior[i];
        }
    }
    
    // Get initial condition name
    virtual std::string name() const = 0;
    
    // Get initial condition type ID
    virtual int type_id() const = 0;
    
    // Evaluate function at a point
    virtual Real evaluate(Real x) const = 0;
};

// ===================================================================
// Step function initial condition
// ===================================================================
class StepInitialCondition : public InitialCondition {
private:
    Real low_value_;
    Real high_value_;
    Real step_start_;
    Real step_end_;
    
public:
    StepInitialCondition(Real low_val = 1.0, Real high_val = 2.0,
                        Real start = 0.5, Real end = 1.0)
        : low_value_(low_val), high_value_(high_val),
          step_start_(start), step_end_(end) {}
    
    ~StepInitialCondition() override = default;
    
    void initialize(Vector& u_interior, const ComputationalDomain& domain) const override {
        const auto& xcc = domain.mesh().cell_centers();
        
        for (int i = 0; i < domain.interior_cells(); ++i) {
            if (xcc[i] >= step_start_ && xcc[i] <= step_end_) {
                u_interior[i] = high_value_;
            } else {
                u_interior[i] = low_value_;
            }
        }
    }
    
    Real evaluate(Real x) const override {
        if (x >= step_start_ && x <= step_end_) {
            return high_value_;
        } else {
            return low_value_;
        }
    }
    
    std::string name() const override {
        return "Step Function (low=" + std::to_string(low_value_) + 
               ", high=" + std::to_string(high_value_) + 
               ", start=" + std::to_string(step_start_) + 
               ", end=" + std::to_string(step_end_) + ")";
    }
    
    int type_id() const override { return 0; }
    
    Real low_value() const { return low_value_; }
    Real high_value() const { return high_value_; }
    Real step_start() const { return step_start_; }
    Real step_end() const { return step_end_; }
};

// ===================================================================
// Sinusoidal initial condition
// ===================================================================
class SinusoidalInitialCondition : public InitialCondition {
private:
    Real amplitude_;
    Real frequency_;
    Real mean_value_;
    Real phase_;
    
public:
    SinusoidalInitialCondition(Real amplitude = 0.5, Real frequency = 1.0,
                              Real mean = 1.0, Real phase = 0.0)
        : amplitude_(amplitude), frequency_(frequency),
          mean_value_(mean), phase_(phase) {}
    
    ~SinusoidalInitialCondition() override = default;
    
    void initialize(Vector& u_interior, const ComputationalDomain& domain) const override {
        const auto& xcc = domain.mesh().cell_centers();
        
        for (int i = 0; i < domain.interior_cells(); ++i) {
            u_interior[i] = mean_value_ + amplitude_ * 
                           std::sin(2.0 * M_PI * frequency_ * xcc[i] + phase_);
        }
    }
    
    Real evaluate(Real x) const override {
        return mean_value_ + amplitude_ * std::sin(2.0 * M_PI * frequency_ * x + phase_);
    }
    
    std::string name() const override {
        return "Sinusoidal (A=" + std::to_string(amplitude_) + 
               ", f=" + std::to_string(frequency_) + 
               ", mean=" + std::to_string(mean_value_) + 
               ", phase=" + std::to_string(phase_) + ")";
    }
    
    int type_id() const override { return 1; }
    
    Real amplitude() const { return amplitude_; }
    Real frequency() const { return frequency_; }
    Real mean_value() const { return mean_value_; }
    Real phase() const { return phase_; }
};

// ===================================================================
// Gaussian initial condition
// ===================================================================
class GaussianInitialCondition : public InitialCondition {
private:
    Real amplitude_;
    Real center_;
    Real width_;
    Real background_;
    
public:
    GaussianInitialCondition(Real amplitude = 1.0, Real center = 1.0,
                            Real width = 0.1, Real background = 1.0)
        : amplitude_(amplitude), center_(center),
          width_(width), background_(background) {}
    
    ~GaussianInitialCondition() override = default;
    
    void initialize(Vector& u_interior, const ComputationalDomain& domain) const override {
        const auto& xcc = domain.mesh().cell_centers();
        
        for (int i = 0; i < domain.interior_cells(); ++i) {
            Real dx = xcc[i] - center_;
            u_interior[i] = background_ + amplitude_ * 
                           std::exp(-dx * dx / (2.0 * width_ * width_));
        }
    }
    
    Real evaluate(Real x) const override {
        Real dx = x - center_;
        return background_ + amplitude_ * std::exp(-dx * dx / (2.0 * width_ * width_));
    }
    
    std::string name() const override {
        return "Gaussian (A=" + std::to_string(amplitude_) + 
               ", center=" + std::to_string(center_) + 
               ", width=" + std::to_string(width_) + 
               ", background=" + std::to_string(background_) + ")";
    }
    
    int type_id() const override { return 2; }
    
    Real amplitude() const { return amplitude_; }
    Real center() const { return center_; }
    Real width() const { return width_; }
    Real background() const { return background_; }
};

// ===================================================================
// Custom function initial condition
// ===================================================================
class CustomInitialCondition : public InitialCondition {
private:
    std::function<Real(Real)> init_func_;
    std::string name_;
    
public:
    CustomInitialCondition(const std::function<Real(Real)>& func,
                          const std::string& name = "Custom")
        : init_func_(func), name_(name) {}
    
    ~CustomInitialCondition() override = default;
    
    void initialize(Vector& u_interior, const ComputationalDomain& domain) const override {
        const auto& xcc = domain.mesh().cell_centers();
        
        for (int i = 0; i < domain.interior_cells(); ++i) {
            u_interior[i] = init_func_(xcc[i]);
        }
    }
    
    Real evaluate(Real x) const override {
        return init_func_(x);
    }
    
    std::string name() const override { return name_; }
    int type_id() const override { return 3; }
};

// ===================================================================
// Initial condition factory
// ===================================================================
class InitialConditionFactory {
public:
    static std::unique_ptr<InitialCondition> create_initial_condition(
        const std::string& ic_type,
        const std::vector<Real>& params = {}) {
        
        if (ic_type == "step" || ic_type == "Step") {
            if (params.size() >= 4) {
                return std::make_unique<StepInitialCondition>(
                    params[0], params[1], params[2], params[3]);
            } else if (params.size() >= 2) {
                return std::make_unique<StepInitialCondition>(
                    params[0], params[1]);
            } else {
                return std::make_unique<StepInitialCondition>();
            }
        }
        else if (ic_type == "sinusoidal" || ic_type == "Sinusoidal") {
            if (params.size() >= 4) {
                return std::make_unique<SinusoidalInitialCondition>(
                    params[0], params[1], params[2], params[3]);
            } else if (params.size() >= 3) {
                return std::make_unique<SinusoidalInitialCondition>(
                    params[0], params[1], params[2]);
            } else if (params.size() >= 2) {
                return std::make_unique<SinusoidalInitialCondition>(
                    params[0], params[1]);
            } else {
                return std::make_unique<SinusoidalInitialCondition>();
            }
        }
        else if (ic_type == "gaussian" || ic_type == "Gaussian") {
            if (params.size() >= 4) {
                return std::make_unique<GaussianInitialCondition>(
                    params[0], params[1], params[2], params[3]);
            } else if (params.size() >= 3) {
                return std::make_unique<GaussianInitialCondition>(
                    params[0], params[1], params[2]);
            } else if (params.size() >= 2) {
                return std::make_unique<GaussianInitialCondition>(
                    params[0], params[1]);
            } else {
                return std::make_unique<GaussianInitialCondition>();
            }
        }
        else {
            throw std::invalid_argument("Unknown initial condition: " + ic_type);
        }
    }
    
    static std::unique_ptr<InitialCondition> create_from_config(
        const CfdConfig& config) {
        // Default to step function for 1D convection
        return std::make_unique<StepInitialCondition>();
    }
    
    static std::vector<std::string> available_initial_conditions() {
        return {
            "Step Function",
            "Sinusoidal",
            "Gaussian",
            "Custom"
        };
    }
};

} // namespace cfd

#endif // INITIALIZER_HPP