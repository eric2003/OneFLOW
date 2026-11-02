#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "domain.hpp"
#include <vector>
#include <memory>
#include <string>
#include <functional>

namespace cfd {

// Type aliases
using Real = double;
using Vector = std::vector<Real>;

// ===================================================================
// Boundary condition base class
// ===================================================================
class BoundaryCondition {
public:
    virtual ~BoundaryCondition() = default;
    
    // Apply boundary condition to solution field
    virtual void apply(Vector& u, const ComputationalDomain& domain) const = 0;
    
    // Apply boundary condition to residual field (if needed)
    virtual void apply_to_residual(Vector& res, const ComputationalDomain& domain) const {
        // Default implementation does nothing
    }
    
    // Get boundary condition type name
    virtual std::string name() const = 0;
    
    // Get boundary condition type ID
    virtual int type_id() const = 0;
    
    // Check if boundary condition is periodic
    virtual bool is_periodic() const { return false; }
    
    // Get ghost cell update method
    virtual void update_ghost_cells(Vector& u, const ComputationalDomain& domain) const = 0;
};

// ===================================================================
// Periodic boundary condition
// ===================================================================
class PeriodicBoundary : public BoundaryCondition {
public:
    PeriodicBoundary() = default;
    ~PeriodicBoundary() override = default;
    
    void apply(Vector& u, const ComputationalDomain& domain) const override {
        domain.apply_periodic_boundary(u);
    }
    
    void update_ghost_cells(Vector& u, const ComputationalDomain& domain) const override {
        // 对于周期边界，直接调用 apply
        apply(u, domain);
    }	
    
    std::string name() const override { return "Periodic"; }
    int type_id() const override { return 0; }
    bool is_periodic() const override { return true; }
};

// ===================================================================
// Dirichlet boundary condition (fixed value)
// ===================================================================
class DirichletBoundary : public BoundaryCondition {
private:
    Real left_value_;
    Real right_value_;
    
public:
    DirichletBoundary(Real left_val = 0.0, Real right_val = 0.0)
        : left_value_(left_val), right_value_(right_val) {}
    
    ~DirichletBoundary() override = default;
    
    void apply(Vector& u, const ComputationalDomain& domain) const override {
        const int nghosts = domain.nghosts();
        const int ist = domain.ist();
        const int ied = domain.ied();
        
        // Left boundary
        for (int i = 0; i < nghosts; ++i) {
            u[i] = left_value_;
        }
        
        // Right boundary
        for (int i = ied + 1; i < domain.ntcells(); ++i) {
            u[i] = right_value_;
        }
    }
    
    void update_ghost_cells(Vector& u, const ComputationalDomain& domain) const override {
        apply(u, domain);
    }
    
    std::string name() const override { 
        return "Dirichlet (left=" + std::to_string(left_value_) + 
               ", right=" + std::to_string(right_value_) + ")";
    }
    
    int type_id() const override { return 1; }
    
    Real left_value() const { return left_value_; }
    Real right_value() const { return right_value_; }
};

// ===================================================================
// Neumann boundary condition (zero gradient)
// ===================================================================
class NeumannBoundary : public BoundaryCondition {
public:
    NeumannBoundary() = default;
    ~NeumannBoundary() override = default;
    
    void apply(Vector& u, const ComputationalDomain& domain) const override {
        const int nghosts = domain.nghosts();
        const int ist = domain.ist();
        const int ied = domain.ied();
        
        // Left boundary: u_ghost = u_interior (zero gradient)
        for (int i = 0; i < nghosts; ++i) {
            u[i] = u[ist + i];
        }
        
        // Right boundary: u_ghost = u_interior (zero gradient)
        for (int i = 0; i < nghosts; ++i) {
            u[ied + 1 + i] = u[ied - nghosts + 1 + i];
        }
    }
    
    void update_ghost_cells(Vector& u, const ComputationalDomain& domain) const override {
        apply(u, domain);
    }
    
    std::string name() const override { return "Neumann (zero gradient)"; }
    int type_id() const override { return 2; }
};

// ===================================================================
// Outflow boundary condition (extrapolation)
// ===================================================================
class OutflowBoundary : public BoundaryCondition {
public:
    OutflowBoundary() = default;
    ~OutflowBoundary() override = default;
    
    void apply(Vector& u, const ComputationalDomain& domain) const override {
        const int nghosts = domain.nghosts();
        const int ist = domain.ist();
        const int ied = domain.ied();
        
        // Left boundary: linear extrapolation
        if (nghosts >= 1) {
            Real slope = u[ist] - u[ist + 1];
            for (int i = nghosts - 1; i >= 0; --i) {
                int ghost_idx = nghosts - 1 - i;
                u[ghost_idx] = u[ist] + (nghosts - i) * slope;
            }
        }
        
        // Right boundary: linear extrapolation
        if (nghosts >= 1) {
            Real slope = u[ied] - u[ied - 1];
            for (int i = 0; i < nghosts; ++i) {
                int ghost_idx = ied + 1 + i;
                u[ghost_idx] = u[ied] + (i + 1) * slope;
            }
        }
    }
    
    void update_ghost_cells(Vector& u, const ComputationalDomain& domain) const override {
        apply(u, domain);
    }
    
    std::string name() const override { return "Outflow (extrapolation)"; }
    int type_id() const override { return 3; }
};

// ===================================================================
// Boundary condition factory
// ===================================================================
class BoundaryFactory {
public:
    static std::unique_ptr<BoundaryCondition> create_boundary(
        const std::string& bc_type, 
        Real left_val = 0.0, 
        Real right_val = 0.0) {
        
        if (bc_type == "periodic" || bc_type == "Periodic") {
            return std::make_unique<PeriodicBoundary>();
        }
        else if (bc_type == "dirichlet" || bc_type == "Dirichlet") {
            return std::make_unique<DirichletBoundary>(left_val, right_val);
        }
        else if (bc_type == "neumann" || bc_type == "Neumann") {
            return std::make_unique<NeumannBoundary>();
        }
        else if (bc_type == "outflow" || bc_type == "Outflow") {
            return std::make_unique<OutflowBoundary>();
        }
        else {
            throw std::invalid_argument("Unknown boundary condition: " + bc_type);
        }
    }
    
    static std::unique_ptr<BoundaryCondition> create_boundary(
        const CfdConfig& config) {
        // Default to periodic for 1D convection
        return std::make_unique<PeriodicBoundary>();
    }
    
    static std::vector<std::string> available_boundary_conditions() {
        return {
            "Periodic",
            "Dirichlet",
            "Neumann",
            "Outflow"
        };
    }
};

} // namespace cfd

#endif // BOUNDARY_HPP