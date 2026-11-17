#include "initial_condition/GaussianInitialCondition.h"
#include <cmath>
#include <sstream>

namespace cfd {

GaussianInitialCondition::GaussianInitialCondition(Real amplitude, 
                                                  Real center,
                                                  Real width, 
                                                  Real background)
    : amplitude_(amplitude), center_(center),
      width_(width), background_(background) {}

void GaussianInitialCondition::initialize(Vector& u_interior, 
                                         const ComputationalDomain& domain) const {
    const auto& xcc = domain.mesh().cell_centers();
    
    for (int i = 0; i < domain.interior_cells(); ++i) {
        Real dx = xcc[i] - center_;
        u_interior[i] = background_ + amplitude_ * 
                       std::exp(-dx * dx / (2.0 * width_ * width_));
    }
}

Real GaussianInitialCondition::evaluate(Real x) const {
    Real dx = x - center_;
    return background_ + amplitude_ * std::exp(-dx * dx / (2.0 * width_ * width_));
}

std::string GaussianInitialCondition::name() const {
    std::ostringstream oss;
    oss << "Gaussian (A=" << amplitude_ 
        << ", center=" << center_ 
        << ", width=" << width_ 
        << ", background=" << background_ << ")";
    return oss.str();
}

std::unique_ptr<InitialCondition> GaussianInitialCondition::clone() const {
    return std::make_unique<GaussianInitialCondition>(*this);
}

} // namespace cfd