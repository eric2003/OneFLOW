#ifndef STEP_INITIAL_CONDITION_H
#define STEP_INITIAL_CONDITION_H

#include "InitialCondition.h"

namespace cfd {

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
                        Real start = 0.5, Real end = 1.0);
    
    ~StepInitialCondition() override = default;
    
    void initialize(Vector& u_interior, const ComputationalDomain& domain) const override;
    
    Real evaluate(Real x) const override;
    
    std::string name() const override;
    
    int type_id() const override { return 0; }
    
    std::unique_ptr<InitialCondition> clone() const override;
    
    // Getters
    Real low_value() const { return low_value_; }
    Real high_value() const { return high_value_; }
    Real step_start() const { return step_start_; }
    Real step_end() const { return step_end_; }
    
    // Setters
    void set_low_value(Real val) { low_value_ = val; }
    void set_high_value(Real val) { high_value_ = val; }
    void set_step_start(Real val) { step_start_ = val; }
    void set_step_end(Real val) { step_end_ = val; }
};

} // namespace cfd

#endif // STEP_INITIAL_CONDITION_H