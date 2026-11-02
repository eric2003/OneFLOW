#ifndef INITIAL_CONDITION_H
#define INITIAL_CONDITION_H

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
    virtual void initialize_with_ghosts(Vector& u_with_ghosts, const ComputationalDomain& domain) const;
    
    // Get initial condition name
    virtual std::string name() const = 0;
    
    // Get initial condition type ID
    virtual int type_id() const = 0;
    
    // Evaluate function at a point
    virtual Real evaluate(Real x) const = 0;
    
    // Clone method for polymorphic copying
    virtual std::unique_ptr<InitialCondition> clone() const = 0;
};

} // namespace cfd

#endif // INITIAL_CONDITION_H