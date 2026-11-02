// ==================== include/mesh.hpp ====================
#ifndef MESH_HPP
#define MESH_HPP

#include "config.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <stdexcept>

namespace cfd {

// Type aliases
using Real = double;
using Vector = std::vector<Real>;

// ===================================================================
// Mesh class for 1D structured grid
// ===================================================================
class Mesh {
private:
    Real xmin_ = 0.0;
    Real xmax_ = 2.0;
    int ncells_ = 40;
    int nnodes_ = 0;
    Real L_ = 0.0;
    Real dx_ = 0.0;
    
    Vector x_;     // Node coordinates
    Vector xcc_;   // Cell center coordinates
    
public:
    // Constructors
    Mesh() = default;
    
    // Constructor with configuration
    explicit Mesh(const CfdConfig& config) 
        : xmin_(config.xmin), xmax_(config.xmax), ncells_(config.ncells) {
        initialize();
    }
    
    // Constructor with explicit parameters
    Mesh(Real xmin, Real xmax, int ncells) 
        : xmin_(xmin), xmax_(xmax), ncells_(ncells) {
        initialize();
    }
    
    // Initialization
    void initialize() {
        nnodes_ = ncells_ + 1;
        L_ = xmax_ - xmin_;
        dx_ = L_ / static_cast<Real>(ncells_);
        
        // Allocate memory
        x_.resize(nnodes_);
        xcc_.resize(ncells_);
        
        // Generate node coordinates
        for (int i = 0; i < nnodes_; ++i) {
            x_[i] = xmin_ + static_cast<Real>(i) * dx_;
        }
        
        // Generate cell center coordinates
        for (int i = 0; i < ncells_; ++i) {
            xcc_[i] = 0.5 * (x_[i] + x_[i + 1]);
        }
    }
    
    // Getters
    Real xmin() const { return xmin_; }
    Real xmax() const { return xmax_; }
    int ncells() const { return ncells_; }
    int nnodes() const { return nnodes_; }
    Real L() const { return L_; }
    Real dx() const { return dx_; }
    
    const Vector& nodes() const { return x_; }
    const Vector& cell_centers() const { return xcc_; }
    
    Real node_coordinate(int i) const {
        if (i < 0 || i >= nnodes_) {
            throw std::out_of_range("Node index out of range");
        }
        return x_[i];
    }
    
    Real cell_center_coordinate(int i) const {
        if (i < 0 || i >= ncells_) {
            throw std::out_of_range("Cell index out of range");
        }
        return xcc_[i];
    }
    
    // Utility methods
    
    // Find cell index containing point x
    int find_cell_index(Real x) const {
        if (x < xmin_ || x > xmax_) {
            throw std::out_of_range("Point outside mesh domain");
        }
        
        int idx = static_cast<int>((x - xmin_) / dx_);
        
        // Handle boundary cases
        if (idx >= ncells_) idx = ncells_ - 1;
        if (idx < 0) idx = 0;
        
        return idx;
    }
    
    // Interpolate value at point x using cell-centered values
    Real interpolate(Real x, const Vector& cell_values) const {
        if (cell_values.size() != static_cast<size_t>(ncells_)) {
            throw std::invalid_argument("Cell values size does not match mesh");
        }
        
        int idx = find_cell_index(x);
        Real x_left = xcc_[idx] - 0.5 * dx_;
        Real x_right = xcc_[idx] + 0.5 * dx_;
        
        // Linear interpolation within the cell
        if (idx < ncells_ - 1) {
            Real alpha = (x - x_left) / dx_;
            return (1.0 - alpha) * cell_values[idx] + alpha * cell_values[idx + 1];
        } else {
            // For last cell, use nearest neighbor
            return cell_values[idx];
        }
    }
    
    // Output and visualization
    void print_info() const {
        std::cout << "Mesh Information:" << std::endl;
        std::cout << "  Domain: [" << xmin_ << ", " << xmax_ << "]" << std::endl;
        std::cout << "  Cells: " << ncells_ << std::endl;
        std::cout << "  Nodes: " << nnodes_ << std::endl;
        std::cout << "  Cell size (dx): " << dx_ << std::endl;
        std::cout << "  Domain length: " << L_ << std::endl;
    }
    
    void print_detailed_info() const {
        print_info();
        std::cout << "\nFirst 5 cells:" << std::endl;
        std::cout << "  Index | Cell Center | Left Node | Right Node" << std::endl;
        std::cout << "  ---------------------------------------------" << std::endl;
        
        for (int i = 0; i < std::min(5, ncells_); ++i) {
            std::cout << "  " << std::setw(5) << i 
                      << " | " << std::setw(11) << std::setprecision(6) << xcc_[i]
                      << " | " << std::setw(9) << std::setprecision(6) << x_[i]
                      << " | " << std::setw(10) << std::setprecision(6) << x_[i + 1] 
                      << std::endl;
        }
        
        if (ncells_ > 5) {
            std::cout << "  ... (showing first 5 of " << ncells_ << " cells)" << std::endl;
        }
    }
    
    // Write mesh data to file
    void write_to_file(const std::string& filename) const {
        std::ofstream file(filename);
        file << "# Mesh Data" << std::endl;
        file << "# xmin = " << xmin_ << std::endl;
        file << "# xmax = " << xmax_ << std::endl;
        file << "# ncells = " << ncells_ << std::endl;
        file << "# dx = " << dx_ << std::endl;
        file << "#" << std::endl;
        file << "# Index, Node_X, CellCenter_X" << std::endl;
        
        file << std::fixed << std::setprecision(6);
        for (int i = 0; i < ncells_; ++i) {
            file << std::setw(6) << i 
                 << std::setw(12) << x_[i]
                 << std::setw(12) << xcc_[i] 
                 << std::endl;
        }
        
        // Last node
        file << std::setw(6) << ncells_
             << std::setw(12) << x_[ncells_]
             << std::setw(12) << "N/A" 
             << std::endl;
        
        file.close();
    }
    
    // Factory methods
    static Mesh create_uniform_mesh(Real xmin, Real xmax, int ncells) {
        return Mesh(xmin, xmax, ncells);
    }
    
    static Mesh create_from_config(const CfdConfig& config) {
        return Mesh(config);
    }
    
    // Create refined mesh (n times more cells)
    Mesh create_refined_mesh(int refinement_factor) const {
        if (refinement_factor <= 0) {
            throw std::invalid_argument("Refinement factor must be positive");
        }
        return Mesh(xmin_, xmax_, ncells_ * refinement_factor);
    }
};

} // namespace cfd

#endif // MESH_HPP