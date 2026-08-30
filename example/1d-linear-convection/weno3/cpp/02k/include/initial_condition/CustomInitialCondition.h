#ifndef CUSTOM_INITIAL_CONDITION_H
#define CUSTOM_INITIAL_CONDITION_H

#include "InitialCondition.h"
#include <functional>

namespace cfd {

// ===================================================================
// Custom function initial condition
// ===================================================================
class CustomInitialCondition : public InitialCondition {
private:
    std::function<Real(Real)> init_func_;
    std::string name_;
    
public:
    CustomInitialCondition(const std::function<Real(Real)>& func,
                          const std::string& name = "Custom");
    
    ~CustomInitialCondition() override = default;
    
    void initialize(Vector& u_interior, const ComputationalDomain& domain) const override;
    
    Real evaluate(Real x) const override;
    
    std::string name() const override;
    
    int type_id() const override { return 3; }
    
    std::unique_ptr<InitialCondition> clone() const override;
    
    // Setter
    void set_name(const std::string& name) { name_ = name; }
};

} // namespace cfd

#endif // CUSTOM_INITIAL_CONDITION_H