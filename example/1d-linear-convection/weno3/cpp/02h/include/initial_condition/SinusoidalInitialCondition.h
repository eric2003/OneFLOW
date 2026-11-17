#ifndef SINUSOIDAL_INITIAL_CONDITION_H
#define SINUSOIDAL_INITIAL_CONDITION_H

#include "InitialCondition.h"

namespace cfd {

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
                              Real mean = 1.0, Real phase = 0.0);
    
    ~SinusoidalInitialCondition() override = default;
    
    void initialize(Vector& u_interior, const ComputationalDomain& domain) const override;
    
    Real evaluate(Real x) const override;
    
    std::string name() const override;
    
    int type_id() const override { return 1; }
    
    std::unique_ptr<InitialCondition> clone() const override;
    
    // Getters
    Real amplitude() const { return amplitude_; }
    Real frequency() const { return frequency_; }
    Real mean_value() const { return mean_value_; }
    Real phase() const { return phase_; }
    
    // Setters
    void set_amplitude(Real val) { amplitude_ = val; }
    void set_frequency(Real val) { frequency_ = val; }
    void set_mean_value(Real val) { mean_value_ = val; }
    void set_phase(Real val) { phase_ = val; }
};

} // namespace cfd

#endif // SINUSOIDAL_INITIAL_CONDITION_H