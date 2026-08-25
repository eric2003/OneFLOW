// ==================== include/flux.hpp ====================
#ifndef FLUX_HPP
#define FLUX_HPP

#include "config.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <memory>
#include <string>
#include <iomanip>
#include <iostream>

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
};

// ===================================================================
// Rusanov (Local Lax-Friedrichs) flux
// ===================================================================
class RusanovFlux : public FluxCalculator {
public:
    RusanovFlux() = default;
    ~RusanovFlux() override = default;
    
    Real compute_flux(Real uL, Real uR, Real wave_speed) const override;
    void compute_flux_array(const Vector& uL, const Vector& uR, 
                           Vector& flux, Real wave_speed) const override;
    
    std::string name() const override { return "Rusanov (Local Lax-Friedrichs)"; }
    int type_id() const override { return 0; }
};

// ===================================================================
// Engquist-Osher flux
// ===================================================================
class EngquistOsherFlux : public FluxCalculator {
public:
    EngquistOsherFlux() = default;
    ~EngquistOsherFlux() override = default;
    
    Real compute_flux(Real uL, Real uR, Real wave_speed) const override;
    void compute_flux_array(const Vector& uL, const Vector& uR, 
                           Vector& flux, Real wave_speed) const override;
    
    std::string name() const override { return "Engquist-Osher"; }
    int type_id() const override { return 1; }
};

// ===================================================================
// Lax-Wendroff flux (for comparison)
// ===================================================================
class LaxWendroffFlux : public FluxCalculator {
public:
    LaxWendroffFlux() = default;
    ~LaxWendroffFlux() override = default;
    
    Real compute_flux(Real uL, Real uR, Real wave_speed) const override;
    void compute_flux_array(const Vector& uL, const Vector& uR, 
                           Vector& flux, Real wave_speed) const override;
    
    std::string name() const override { return "Lax-Wendroff"; }
    int type_id() const override { return 2; }
};

// ===================================================================
// Lax-Wendroff flux with dt and dx parameters
// ===================================================================
class LaxWendroffFluxExtended : public LaxWendroffFlux {
private:
    Real dt_ = 0.025;
    Real dx_ = 1.0;
    
public:
    LaxWendroffFluxExtended(Real dt, Real dx) : dt_(dt), dx_(dx) {}
    
    Real compute_flux(Real uL, Real uR, Real wave_speed) const override;
    void compute_flux_array(const Vector& uL, const Vector& uR, 
                           Vector& flux, Real wave_speed) const override;
    
    std::string name() const override;
    int type_id() const override { return 2; }
};

// ===================================================================
// Upwind flux (first-order)
// ===================================================================
class UpwindFlux : public FluxCalculator {
public:
    UpwindFlux() = default;
    ~UpwindFlux() override = default;
    
    Real compute_flux(Real uL, Real uR, Real wave_speed) const override;
    void compute_flux_array(const Vector& uL, const Vector& uR, 
                           Vector& flux, Real wave_speed) const override;
    
    std::string name() const override { return "Upwind"; }
    int type_id() const override { return 3; }
};

// ===================================================================
// Central flux (unstable for hyperbolic equations, but included for reference)
// ===================================================================
class CentralFlux : public FluxCalculator {
public:
    CentralFlux() = default;
    ~CentralFlux() override = default;
    
    Real compute_flux(Real uL, Real uR, Real wave_speed) const override;
    void compute_flux_array(const Vector& uL, const Vector& uR, 
                           Vector& flux, Real wave_speed) const override;
    
    std::string name() const override { return "Central"; }
    int type_id() const override { return 4; }
};

// ===================================================================
// Roe flux
// ===================================================================
class RoeFlux : public FluxCalculator {
public:
    RoeFlux() = default;
    ~RoeFlux() override = default;
    
    Real compute_flux(Real uL, Real uR, Real wave_speed) const override;
    void compute_flux_array(const Vector& uL, const Vector& uR, 
                           Vector& flux, Real wave_speed) const override;
    
    std::string name() const override { return "Roe"; }
    int type_id() const override { return 5; }
};

// ===================================================================
// HLL flux
// ===================================================================
class HLLFlux : public FluxCalculator {
public:
    HLLFlux() = default;
    ~HLLFlux() override = default;
    
    Real compute_flux(Real uL, Real uR, Real wave_speed) const override;
    void compute_flux_array(const Vector& uL, const Vector& uR, 
                           Vector& flux, Real wave_speed) const override;
    
    std::string name() const override { return "HLL"; }
    int type_id() const override { return 6; }
};

// ===================================================================
// Flux factory
// ===================================================================
class FluxFactory {
public:
    // Basic factory methods
    static std::unique_ptr<FluxCalculator> create_flux_calculator(int flux_type);
    static std::unique_ptr<FluxCalculator> create_flux_calculator(const std::string& flux_name);
    
    // Extended factory methods with dt and dx parameters
    static std::unique_ptr<FluxCalculator> create_flux_calculator(int flux_type, Real dt, Real dx);
    static std::unique_ptr<FluxCalculator> create_flux_calculator(const std::string& flux_name, Real dt, Real dx);
    
    // Extended factory with all flux types
    static std::unique_ptr<FluxCalculator> create_extended_flux_calculator(int flux_type, Real dt, Real dx);
    
    // Get available fluxes
    static std::vector<std::string> available_fluxes();
    static std::vector<std::string> available_fluxes_extended();
    
    // Flux comparison utilities
    static Vector compute_flux_comparison(Real uL, Real uR, Real wave_speed);
    static void print_flux_comparison(Real uL, Real uR, Real wave_speed);
};

} // namespace cfd

#endif // FLUX_HPP