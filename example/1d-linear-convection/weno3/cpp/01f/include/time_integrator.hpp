#ifndef TIME_INTEGRATOR_HPP
#define TIME_INTEGRATOR_HPP

#include "residual.hpp"
#include "domain.hpp"
#include <memory>

namespace cfd {

// Type aliases
using Real = double;
using Vector = std::vector<Real>;

// ===================================================================
// Time integrator base class
// ===================================================================
class TimeIntegrator {
protected:
    std::unique_ptr<ResidualCalculator> residual_calculator_;
    
public:
    TimeIntegrator(std::unique_ptr<ResidualCalculator> residual_calculator)
        : residual_calculator_(std::move(residual_calculator)) {}
    
    virtual ~TimeIntegrator() = default;
    
    // Perform one time step
    virtual void step(Vector& u_with_ghosts, 
                     Real dt,
                     const ComputationalDomain& domain) const = 0;
    
    // Get time integrator name
    virtual std::string name() const = 0;
    
    // Get order of accuracy
    virtual int order() const = 0;
    
    // Get residual calculator
    const ResidualCalculator& residual_calculator() const { 
        return *residual_calculator_; 
    }
};

// ===================================================================
// Forward Euler (RK1) time integrator
// ===================================================================
class ForwardEulerIntegrator : public TimeIntegrator {
public:
    ForwardEulerIntegrator(std::unique_ptr<ResidualCalculator> residual_calculator)
        : TimeIntegrator(std::move(residual_calculator)) {}
    
    ~ForwardEulerIntegrator() override = default;
    
    void step(Vector& u_with_ghosts, 
             Real dt,
             const ComputationalDomain& domain) const override {
        
        // Get domain parameters
        const int ncells = domain.mesh().ncells();
        const int ist = domain.ist();
        Real dx = domain.mesh().dx();
        
        // Compute residual
        Vector residual(ncells);
        residual_calculator_->compute(u_with_ghosts, residual, domain, dx);
        
        // Update solution (only interior cells)
        for (int i = 0; i < ncells; ++i) {
            u_with_ghosts[ist + i] += dt * residual[i];
        }
    }
    
    std::string name() const override { 
        return "Forward Euler (RK1)";
    }
    
    int order() const override { return 1; }
};

// ===================================================================
// RK2 time integrator (midpoint method)
// ===================================================================
class RK2Integrator : public TimeIntegrator {
public:
    RK2Integrator(std::unique_ptr<ResidualCalculator> residual_calculator)
        : TimeIntegrator(std::move(residual_calculator)) {}
    
    ~RK2Integrator() override = default;
    
    void step(Vector& u_with_ghosts, 
             Real dt,
             const ComputationalDomain& domain) const override {
        
        // Get domain parameters
        const int ncells = domain.mesh().ncells();
        const int ist = domain.ist();
        Real dx = domain.mesh().dx();
        
        // Stage 1
        Vector residual(ncells);
        residual_calculator_->compute(u_with_ghosts, residual, domain, dx);
        
        // Temporary storage for stage 1 solution
        Vector u_stage = u_with_ghosts;
        
        // Update to stage 1
        for (int i = 0; i < ncells; ++i) {
            u_stage[ist + i] = u_with_ghosts[ist + i] + dt * residual[i];
        }
        
        // Stage 2
        residual_calculator_->compute(u_stage, residual, domain, dx);
        
        // Final update
        for (int i = 0; i < ncells; ++i) {
            u_with_ghosts[ist + i] = 0.5 * u_with_ghosts[ist + i] + 
                                    0.5 * u_stage[ist + i] + 
                                    0.5 * dt * residual[i];
        }
    }
    
    std::string name() const override { 
        return "Runge-Kutta 2 (Midpoint)";
    }
    
    int order() const override { return 2; }
};

// ===================================================================
// TVD RK3 time integrator
// ===================================================================
class TVDRK3Integrator : public TimeIntegrator {
public:
    TVDRK3Integrator(std::unique_ptr<ResidualCalculator> residual_calculator)
        : TimeIntegrator(std::move(residual_calculator)) {}
    
    ~TVDRK3Integrator() override = default;
    
    void step(Vector& u_with_ghosts, 
             Real dt,
             const ComputationalDomain& domain) const override {
        
        // Get domain parameters
        const int ncells = domain.mesh().ncells();
        const int ist = domain.ist();
        Real dx = domain.mesh().dx();
        
        Vector residual(ncells);
        
        // Stage 1
        residual_calculator_->compute(u_with_ghosts, residual, domain, dx);
        Vector u1 = u_with_ghosts;
        for (int i = 0; i < ncells; ++i) {
            u1[ist + i] = u_with_ghosts[ist + i] + dt * residual[i];
        }
        
        // Stage 2
        residual_calculator_->compute(u1, residual, domain, dx);
        Vector u2 = u_with_ghosts;
        for (int i = 0; i < ncells; ++i) {
            u2[ist + i] = 0.75 * u_with_ghosts[ist + i] + 
                         0.25 * u1[ist + i] + 
                         0.25 * dt * residual[i];
        }
        
        // Stage 3
        residual_calculator_->compute(u2, residual, domain, dx);
        for (int i = 0; i < ncells; ++i) {
            u_with_ghosts[ist + i] = (1.0/3.0) * u_with_ghosts[ist + i] + 
                                    (2.0/3.0) * u2[ist + i] + 
                                    (2.0/3.0) * dt * residual[i];
        }
    }
    
    std::string name() const override { 
        return "TVD Runge-Kutta 3";
    }
    
    int order() const override { return 3; }
};

// ===================================================================
// Time integrator factory
// ===================================================================
class TimeIntegratorFactory {
public:
    static std::unique_ptr<TimeIntegrator> create_integrator(
        const std::string& method,
        std::unique_ptr<ResidualCalculator> residual_calculator) {
        
        if (method == "rk1" || method == "euler" || method == "ForwardEuler") {
            return std::make_unique<ForwardEulerIntegrator>(
                std::move(residual_calculator));
        }
        else if (method == "rk2" || method == "RK2" || method == "midpoint") {
            return std::make_unique<RK2Integrator>(
                std::move(residual_calculator));
        }
        else if (method == "rk3" || method == "RK3" || method == "TVDRK3") {
            return std::make_unique<TVDRK3Integrator>(
                std::move(residual_calculator));
        }
        else {
            throw std::invalid_argument("Unknown time integration method: " + method);
        }
    }
    
    static std::unique_ptr<TimeIntegrator> create_integrator(
        int order,
        std::unique_ptr<ResidualCalculator> residual_calculator) {
        
        switch (order) {
            case 1:
                return std::make_unique<ForwardEulerIntegrator>(
                    std::move(residual_calculator));
            case 2:
                return std::make_unique<RK2Integrator>(
                    std::move(residual_calculator));
            case 3:
                return std::make_unique<TVDRK3Integrator>(
                    std::move(residual_calculator));
            default:
                throw std::invalid_argument("Unsupported RK order: " + std::to_string(order));
        }
    }
    
    static std::vector<std::string> available_integrators() {
        return {
            "Forward Euler (RK1)",
            "Runge-Kutta 2 (RK2)",
            "TVD Runge-Kutta 3 (RK3)"
        };
    }
};

} // namespace cfd

#endif // TIME_INTEGRATOR_HPP