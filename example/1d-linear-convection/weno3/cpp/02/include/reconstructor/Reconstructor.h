#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include "domain.hpp"
#include <vector>
#include <memory>
#include <string>

namespace cfd {

// Type aliases
using Real = double;
using Vector = std::vector<Real>;

// ===================================================================
// Reconstructor base class
// ===================================================================
class Reconstructor {
public:
    virtual ~Reconstructor() = default;
    
    // Reconstruct face values from cell averages
    virtual void reconstruct(const Vector& q, 
                           Vector& q_face_left,
                           Vector& q_face_right,
                           const ComputationalDomain& domain) = 0;
    
    // Get reconstructor name
    virtual std::string name() const = 0;
    
    // Get spatial order
    virtual int order() const = 0;
    
    // Get stencil width (number of cells needed on each side)
    virtual int stencil_width() const = 0;
    
    // Check if reconstructor is linear
    virtual bool is_linear() const { return true; }
    
    // Check if reconstructor is nonlinear (e.g., WENO)
    virtual bool is_nonlinear() const { return !is_linear(); }
    
    // Check if reconstructor is TVD (Total Variation Diminishing)
    virtual bool is_tvd() const { return false; }
    
    // Check if reconstructor is monotonicity preserving
    virtual bool is_monotonic() const { return false; }
    
    // Get ghost cells required
    virtual int required_ghost_cells() const { return stencil_width() - 1; }
    
    // Apply limiter (for schemes with limiters)
    virtual void apply_limiter(Vector& slopes) const {}
    
    // Print reconstructor info
    virtual void print_info() const {
        std::cout << name() << " (Order: " << order() 
                  << ", Linear: " << (is_linear() ? "Yes" : "No")
                  << ", TVD: " << (is_tvd() ? "Yes" : "No")
                  << ", Stencil: " << stencil_width() << ")" << std::endl;
    }
};

} // namespace cfd

#endif // RECONSTRUCTOR_H