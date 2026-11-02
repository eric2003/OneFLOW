#include "initial_condition/StepInitialCondition.h"
#include <sstream>

namespace cfd {

StepInitialCondition::StepInitialCondition(Real low_val, Real high_val,
                                         Real start, Real end)
    : low_value_(low_val), high_value_(high_val),
      step_start_(start), step_end_(end) {}

void StepInitialCondition::initialize(Vector& u_interior, 
                                     const ComputationalDomain& domain) const {
    const auto& xcc = domain.mesh().cell_centers();
    
    for (int i = 0; i < domain.interior_cells(); ++i) {
        if (xcc[i] >= step_start_ && xcc[i] <= step_end_) {
            u_interior[i] = high_value_;
        } else {
            u_interior[i] = low_value_;
        }
    }
}

Real StepInitialCondition::evaluate(Real x) const {
    if (x >= step_start_ && x <= step_end_) {
        return high_value_;
    } else {
        return low_value_;
    }
}

std::string StepInitialCondition::name() const {
    std::ostringstream oss;
    oss << "Step Function (low=" << low_value_ 
        << ", high=" << high_value_ 
        << ", start=" << step_start_ 
        << ", end=" << step_end_ << ")";
    return oss.str();
}

std::unique_ptr<InitialCondition> StepInitialCondition::clone() const {
    return std::make_unique<StepInitialCondition>(*this);
}

} // namespace cfd