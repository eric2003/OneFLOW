#ifndef FLUX_CALCULATOR_H
#define FLUX_CALCULATOR_H

#include <vector>
#include <string>
#include <memory>

namespace cfd {

// Type aliases
using Real = double;
using Vector = std::vector<Real>;

// ===================================================================
// Flux calculator base class
// ===================================================================
class FluxCalculator {
public:
    virtual ~FluxCalculator() = default;
    
    // Compute flux from left and right states
    virtual Real compute_flux(Real uL, Real uR, Real wave_speed) const = 0;
    
    // Compute flux for arrays
    virtual void compute_flux_array(const Vector& uL, const Vector& uR, 
                                   Vector& flux, Real wave_speed) const = 0;
    
    // Get flux type name
    virtual std::string name() const = 0;
    
    // Get flux type ID
    virtual int type_id() const = 0;
    
    // Get flux category
    virtual std::string category() const { return "General"; }
    
    // Check if flux is upwind-based
    virtual bool is_upwind() const { return false; }
    
    // Check if flux is centered
    virtual bool is_centered() const { return false; }
    
    // Check if flux has dissipation
    virtual bool has_dissipation() const { return true; }
    
    // Get dissipation coefficient (if applicable)
    virtual Real dissipation_coefficient() const { return 1.0; }
    
    // Check if flux is conservative
    virtual bool is_conservative() const { return true; }
    
    // Check if flux is entropy stable
    virtual bool is_entropy_stable() const { return false; }
};

} // namespace cfd

#endif // FLUX_CALCULATOR_H