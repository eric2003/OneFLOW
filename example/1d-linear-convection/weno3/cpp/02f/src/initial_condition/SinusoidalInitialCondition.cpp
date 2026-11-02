#include "initial_condition/SinusoidalInitialCondition.h"
#include <cmath>
#include <sstream>

namespace cfd {

SinusoidalInitialCondition::SinusoidalInitialCondition(Real amplitude, 
                                                      Real frequency,
                                                      Real mean, 
                                                      Real phase)
    : amplitude_(amplitude), frequency_(frequency),
      mean_value_(mean), phase_(phase) {}

void SinusoidalInitialCondition::initialize(Vector& u_interior, 
                                           const ComputationalDomain& domain) const {
    const auto& xcc = domain.mesh().cell_centers();
    
    for (int i = 0; i < domain.interior_cells(); ++i) {
        u_interior[i] = mean_value_ + amplitude_ * 
                       std::sin(2.0 * M_PI * frequency_ * xcc[i] + phase_);
    }
}

Real SinusoidalInitialCondition::evaluate(Real x) const {
    return mean_value_ + amplitude_ * std::sin(2.0 * M_PI * frequency_ * x + phase_);
}

std::string SinusoidalInitialCondition::name() const {
    std::ostringstream oss;
    oss << "Sinusoidal (A=" << amplitude_ 
        << ", f=" << frequency_ 
        << ", mean=" << mean_value_ 
        << ", phase=" << phase_ << ")";
    return oss.str();
}

std::unique_ptr<InitialCondition> SinusoidalInitialCondition::clone() const {
    return std::make_unique<SinusoidalInitialCondition>(*this);
}

} // namespace cfd