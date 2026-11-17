#ifndef RESIDUAL_HPP
#define RESIDUAL_HPP

#include "domain.hpp"
#include "solution.hpp"
#include "reconstructor/ReconstructorFactory.h"
#include "reconstructor/advanced/EnoReconstructor.h"
#include "reconstructor/advanced/WenoReconstructor.h"
#include "flux/FluxFactory.h"
#include "flux/basic/RusanovFlux.h"
#include "boundary.hpp"
#include <memory>
#include <vector>

namespace cfd {

// Type aliases
using Real = double;
using Vector = std::vector<Real>;

// ===================================================================
// Residual calculator base class
// ===================================================================
class ResidualCalculator {
protected:
    std::unique_ptr<Reconstructor> reconstructor_;
    std::unique_ptr<FluxCalculator> flux_calculator_;
    std::unique_ptr<BoundaryCondition> boundary_condition_;
    Real wave_speed_;
    
public:
    ResidualCalculator(std::unique_ptr<Reconstructor> reconstructor,
                      std::unique_ptr<FluxCalculator> flux_calculator,
                      std::unique_ptr<BoundaryCondition> boundary_condition,
                      Real wave_speed)
        : reconstructor_(std::move(reconstructor))
        , flux_calculator_(std::move(flux_calculator))
        , boundary_condition_(std::move(boundary_condition))
        , wave_speed_(wave_speed) {}
    
    virtual ~ResidualCalculator() = default;
    
    // Compute residual for the given solution state
    virtual void compute(Vector& u_with_ghosts, 
                        Vector& residual,
                        const ComputationalDomain& domain,
                        Real dx) const = 0;
    
    // Get residual calculator name
    virtual std::string name() const = 0;
    
    // Get components
    const Reconstructor& reconstructor() const { return *reconstructor_; }
    const FluxCalculator& flux_calculator() const { return *flux_calculator_; }
    const BoundaryCondition& boundary_condition() const { return *boundary_condition_; }
    Real wave_speed() const { return wave_speed_; }
};

// ===================================================================
// Standard residual calculator for 1D convection
// ===================================================================
class ConvectionResidualCalculator : public ResidualCalculator {
public:
    ConvectionResidualCalculator(std::unique_ptr<Reconstructor> reconstructor,
                                std::unique_ptr<FluxCalculator> flux_calculator,
                                std::unique_ptr<BoundaryCondition> boundary_condition,
                                Real wave_speed)
        : ResidualCalculator(std::move(reconstructor), 
                           std::move(flux_calculator), 
                           std::move(boundary_condition), 
                           wave_speed) {}
    
    ~ConvectionResidualCalculator() override = default;
    
    void compute(Vector& u_with_ghosts, 
                Vector& residual,
                const ComputationalDomain& domain,
                Real dx) const override {
        
        // Apply boundary conditions
        boundary_condition_->update_ghost_cells(u_with_ghosts, domain);
        
        // Get domain parameters
        const int ncells = domain.mesh().ncells();
        const int ist = domain.ist();
        
        // Temporary arrays for face values
        std::vector<Real> q_face_left(ncells + 1);
        std::vector<Real> q_face_right(ncells + 1);
        std::vector<Real> flux(ncells + 1);
        
        // Reconstruct face values
        reconstructor_->reconstruct(u_with_ghosts, 
                                  q_face_left,
                                  q_face_right,
                                  domain);
        
        // Compute fluxes at faces
        flux_calculator_->compute_flux_array(
            q_face_left,
            q_face_right,
            flux,
            wave_speed_
        );
        
        // Compute residual (flux divergence)
        Real dx_inv = 1.0 / dx;
        for (int i = 0; i < ncells; ++i) {
            residual[i] = -(flux[i+1] - flux[i]) * dx_inv;
        }
    }
    
    std::string name() const override {
        return "Convection Residual (" + 
               reconstructor_->name() + " + " + 
               flux_calculator_->name() + ")";
    }
};

// ===================================================================
// Compact residual calculator (saves memory by reusing arrays)
// ===================================================================
class CompactConvectionResidualCalculator : public ResidualCalculator {
private:
    mutable std::vector<Real> q_face_left_;
    mutable std::vector<Real> q_face_right_;
    mutable std::vector<Real> flux_;
    
public:
    CompactConvectionResidualCalculator(std::unique_ptr<Reconstructor> reconstructor,
                                       std::unique_ptr<FluxCalculator> flux_calculator,
                                       std::unique_ptr<BoundaryCondition> boundary_condition,
                                       Real wave_speed)
        : ResidualCalculator(std::move(reconstructor), 
                           std::move(flux_calculator), 
                           std::move(boundary_condition), 
                           wave_speed) {}
    
    ~CompactConvectionResidualCalculator() override = default;
    
    void compute(Vector& u_with_ghosts, 
                Vector& residual,
                const ComputationalDomain& domain,
                Real dx) const override {
        
        // Apply boundary conditions
        boundary_condition_->update_ghost_cells(u_with_ghosts, domain);
        
        // Get domain parameters
        const int ncells = domain.mesh().ncells();
        
        // Resize temporary arrays if needed
        if (q_face_left_.size() != static_cast<size_t>(ncells + 1)) {
            q_face_left_.resize(ncells + 1);
            q_face_right_.resize(ncells + 1);
            flux_.resize(ncells + 1);
        }
        
        // Reconstruct face values
        reconstructor_->reconstruct(u_with_ghosts, 
                                  q_face_left_,
                                  q_face_right_,
                                  domain);
        
        // Compute fluxes at faces
        flux_calculator_->compute_flux_array(
            q_face_left_,
            q_face_right_,
            flux_,
            wave_speed_
        );
        
        // Compute residual (flux divergence)
        Real dx_inv = 1.0 / dx;
        for (int i = 0; i < ncells; ++i) {
            residual[i] = -(flux_[i+1] - flux_[i]) * dx_inv;
        }
    }
    
    std::string name() const override {
        return "Compact Convection Residual (" + 
               reconstructor_->name() + " + " + 
               flux_calculator_->name() + ")";
    }
};

// ===================================================================
// Residual calculator factory
// ===================================================================
class ResidualCalculatorFactory {
public:
    static std::unique_ptr<ResidualCalculator> create_convection_residual(
        std::unique_ptr<Reconstructor> reconstructor,
        std::unique_ptr<FluxCalculator> flux_calculator,
        std::unique_ptr<BoundaryCondition> boundary_condition,
        Real wave_speed,
        bool compact = false) {
        
        if (compact) {
            return std::make_unique<CompactConvectionResidualCalculator>(
                std::move(reconstructor),
                std::move(flux_calculator),
                std::move(boundary_condition),
                wave_speed
            );
        } else {
            return std::make_unique<ConvectionResidualCalculator>(
                std::move(reconstructor),
                std::move(flux_calculator),
                std::move(boundary_condition),
                wave_speed
            );
        }
    }
    
    static std::unique_ptr<ResidualCalculator> create_from_config(
        const CfdConfig& config,
        bool compact = false) {
        
        // Create reconstructor
        auto reconstructor = ReconstructorFactory::create_reconstructor(config);
        
        // Create flux calculator
        auto flux_calculator = FluxFactory::create_flux_calculator(config.flux_type);
        
        // Create boundary condition (default: periodic)
        auto boundary_condition = BoundaryFactory::create_boundary(config);
        
        return create_convection_residual(
            std::move(reconstructor),
            std::move(flux_calculator),
            std::move(boundary_condition),
            config.wave_speed,
            compact
        );
    }
};

} // namespace cfd

#endif // RESIDUAL_HPP