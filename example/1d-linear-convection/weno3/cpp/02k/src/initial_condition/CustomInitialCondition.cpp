#include "initial_condition/CustomInitialCondition.h"

namespace cfd {

CustomInitialCondition::CustomInitialCondition(const std::function<Real(Real)>& func,
                                              const std::string& name)
    : init_func_(func), name_(name) {}

void CustomInitialCondition::initialize(Vector& u_interior, 
                                       const ComputationalDomain& domain) const {
    const auto& xcc = domain.mesh().cell_centers();
    
    for (int i = 0; i < domain.interior_cells(); ++i) {
        u_interior[i] = init_func_(xcc[i]);
    }
}

Real CustomInitialCondition::evaluate(Real x) const {
    return init_func_(x);
}

std::string CustomInitialCondition::name() const { 
    return name_; 
}

std::unique_ptr<InitialCondition> CustomInitialCondition::clone() const {
    return std::make_unique<CustomInitialCondition>(*this);
}

} // namespace cfd