// ==================== include/domain.hpp ====================
#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include "mesh.hpp"
#include "config.hpp"
#include <vector>
#include <iostream>
#include <stdexcept>
#include <memory>

namespace cfd {

// ===================================================================
// Computational domain class
// ===================================================================
class ComputationalDomain {
private:
    Mesh mesh_;
    CfdConfig config_;
    int nghosts_ = 0;
    int ist_ = 0;     // Start index of interior cells (1-based)
    int ied_ = 0;     // End index of interior cells (1-based)
    int ntcells_ = 0; // Total cells including ghost cells
    
public:
    // ===================================================================
    // Constructors
    // ===================================================================
    ComputationalDomain() = default;
    
    // Constructor with configuration
    explicit ComputationalDomain(const CfdConfig& config) {
        initialize(config);
    }
    
    // Constructor with mesh and configuration
    ComputationalDomain(const Mesh& mesh, const CfdConfig& config) 
        : mesh_(mesh), config_(config) {
        calculate_ghost_cells();
    }
    
    // ===================================================================
    // Initialization
    // ===================================================================
    void initialize(const CfdConfig& config) {
        config_ = config;
        config_.validate();
        
        // Create mesh from configuration
        mesh_ = Mesh(config_);
        
        // Calculate ghost cells
        calculate_ghost_cells();
    }
    
    void calculate_ghost_cells() {
        // Calculate number of ghost cells based on reconstruction scheme
        if (config_.recon_scheme == "eno") {
            nghosts_ = config_.spatial_order;
        } else if (config_.recon_scheme == "weno") {
            nghosts_ = config_.spatial_order / 2 + 1;
        } else {
            throw std::runtime_error("Unknown reconstruction scheme: " + config_.recon_scheme);
        }
        
        // Calculate indices (1-based indexing for compatibility with Fortran)
        ist_ = nghosts_;
        ied_ = ist_ + mesh_.ncells();
        ntcells_ = mesh_.ncells() + 2 * nghosts_;
        
        // Validate indices
        if (ist_ <= 0 || ied_ < ist_) {
            throw std::runtime_error("Invalid domain indices calculated");
        }
    }
    
    // ===================================================================
    // Getters
    // ===================================================================
    const Mesh& mesh() const { return mesh_; }
    const CfdConfig& config() const { return config_; }
    
    int nghosts() const { return nghosts_; }
    int ist() const { return ist_; }
    int ied() const { return ied_; }
    int ntcells() const { return ntcells_; }
    
    // Number of interior cells (physical domain)
    int interior_cells() const { return mesh_.ncells(); }
    
    // ===================================================================
    // Index conversion methods
    // ===================================================================
    
    // Convert global index (with ghosts) to interior index (0-based)
    int global_to_interior(int global_idx) const {
        if (global_idx < ist_ || global_idx > ied_) {
            throw std::out_of_range("Global index is not an interior cell");
        }
        return global_idx - ist_;
    }
    
    // Convert interior index (0-based) to global index (with ghosts)
    int interior_to_global(int interior_idx) const {
        if (interior_idx < 0 || interior_idx >= mesh_.ncells()) {
            throw std::out_of_range("Interior index out of range");
        }
        return ist_ + interior_idx;
    }
    
    // Check if index is in interior domain
    bool is_interior_cell(int global_idx) const {
        return (global_idx >= ist_ && global_idx <= ied_);
    }
    
    // Check if index is ghost cell
    bool is_ghost_cell(int global_idx) const {
        return (global_idx < ist_ || global_idx > ied_);
    }
    
    // ===================================================================
    // Domain operations
    // ===================================================================
    
    // Apply periodic boundary conditions to a field
    template<typename VectorType>
    void apply_periodic_boundary(VectorType& field) const {
        if (field.size() != static_cast<size_t>(ntcells_)) {
            throw std::invalid_argument("Field size does not match domain");
        }
        
        // Left ghost cells = right interior cells
        for (int ig = 0; ig < nghosts_; ++ig) {
            int ghost_idx = ist_ - 1 - ig;
            int interior_idx = ied_ - 1 - ig;
            field[ghost_idx] = field[interior_idx];
        }
        
        // Right ghost cells = left interior cells
        for (int ig = 0; ig < nghosts_; ++ig) {
            int ghost_idx = ied_ + ig;
            int interior_idx = ist_ + ig;
            field[ghost_idx] = field[interior_idx];
        }
    }
    
    // Extract interior solution (without ghost cells)
    template<typename VectorType>
    VectorType extract_interior(const VectorType& field_with_ghosts) const {
        if (field_with_ghosts.size() != static_cast<size_t>(ntcells_)) {
            throw std::invalid_argument("Field size does not match domain");
        }
        
        VectorType interior(mesh_.ncells());
        for (int i = 0; i < mesh_.ncells(); ++i) {
            interior[i] = field_with_ghosts[ist_ + i];
        }
        return interior;
    }
    
    // Create field with ghost cells initialized to zero
    template<typename VectorType>
    VectorType create_field_with_ghosts() const {
        return VectorType(ntcells_, 0.0);
    }
    
    // Create field with ghost cells initialized with interior values
    template<typename VectorType>
    VectorType create_field_with_ghosts(const VectorType& interior_values) const {
        if (interior_values.size() != static_cast<size_t>(mesh_.ncells())) {
            throw std::invalid_argument("Interior values size does not match mesh");
        }
        
        VectorType field = create_field_with_ghosts<VectorType>();
        
        // Copy interior values
        for (int i = 0; i < mesh_.ncells(); ++i) {
            field[ist_ + i] = interior_values[i];
        }
        
        // Apply boundary conditions
        apply_periodic_boundary(field);
        
        return field;
    }
    
    // ===================================================================
    // Validation and diagnostics
    // ===================================================================
    void validate() const {
        if (nghosts_ <= 0) {
            throw std::runtime_error("Invalid number of ghost cells: " + std::to_string(nghosts_));
        }
        
        if (ist_ <= 0 || ied_ < ist_) {
            throw std::runtime_error("Invalid domain indices: ist=" + std::to_string(ist_) + 
                                   ", ied=" + std::to_string(ied_));
        }
        
        if (ntcells_ != mesh_.ncells() + 2 * nghosts_) {
            throw std::runtime_error("Inconsistent total cell count");
        }
    }
    
    void print_info() const {
        std::cout << "Computational Domain:" << std::endl;
        std::cout << "  Interior cells: " << mesh_.ncells() << std::endl;
        std::cout << "  Ghost cells: " << nghosts_ << " (each side)" << std::endl;
        std::cout << "  Total cells: " << ntcells_ << std::endl;
        std::cout << "  Interior range: [" << ist_ << ", " << ied_ << "]" << std::endl;
        std::cout << "  Reconstruction: " << config_.recon_scheme << config_.spatial_order << std::endl;
        std::cout << "  Ghost calculation: ";
        
        if (config_.recon_scheme == "eno") {
            std::cout << "ENO" << config_.spatial_order << " requires " << config_.spatial_order << " ghosts" << std::endl;
        } else if (config_.recon_scheme == "weno") {
            std::cout << "WENO" << config_.spatial_order << " requires " 
                     << (config_.spatial_order / 2 + 1) << " ghosts" << std::endl;
        }
    }
    
    void print_detailed_info() const {
        print_info();
        mesh_.print_info();
        
        std::cout << "\nDomain Layout:" << std::endl;
        std::cout << "  Left Ghost: indices 0 to " << (ist_ - 2) << " (" << (ist_ - 1) << " cells)" << std::endl;
        std::cout << "  Interior: indices " << (ist_ - 1) << " to " << (ied_ - 1) 
                 << " (" << mesh_.ncells() << " cells)" << std::endl;
        std::cout << "  Right Ghost: indices " << ied_ << " to " << (ntcells_ - 1) 
                 << " (" << nghosts_ << " cells)" << std::endl;
    }
    
    // ===================================================================
    // Factory methods
    // ===================================================================
    static ComputationalDomain create_from_config(const CfdConfig& config) {
        return ComputationalDomain(config);
    }
    
    static ComputationalDomain create_with_mesh(const Mesh& mesh, const CfdConfig& config) {
        return ComputationalDomain(mesh, config);
    }
};

} // namespace cfd

#endif // DOMAIN_HPP