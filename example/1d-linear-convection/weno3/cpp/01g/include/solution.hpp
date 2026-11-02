// ==================== include/solution.hpp ====================
#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include "domain.hpp"
#include "initializer.hpp"
#include "boundary.hpp"
#include <vector>
#include <functional>
#include <algorithm>
#include <iomanip>  // for std::setw, std::setprecision
#include <fstream>  // for std::ofstream
#include <numeric>  // for std::accumulate

namespace cfd {

// Type aliases
using Real = double;
using Vector = std::vector<Real>;

// ===================================================================
// Solution class for CFD
// ===================================================================
class Solution {
private:
    ComputationalDomain domain_;
    
    // Solution fields
    Vector q_face_left_;   // Left face values
    Vector q_face_right_;  // Right face values
    Vector flux_;          // Flux at faces
    Vector res_;           // Residual
    Vector u_;             // Current solution with ghost cells
    Vector un_;            // Old solution with ghost cells
    
public:
    // ===================================================================
    // Constructors
    // ===================================================================
    Solution() = default;
    
    explicit Solution(const ComputationalDomain& domain) {
        initialize(domain);
    }
    
    // ===================================================================
    // Initialization
    // ===================================================================
    void initialize(const ComputationalDomain& domain) {
        domain_ = domain;
        
        // Allocate memory
        q_face_left_.resize(domain_.mesh().nnodes(), 0.0);
        q_face_right_.resize(domain_.mesh().nnodes(), 0.0);
        flux_.resize(domain_.mesh().nnodes(), 0.0);
        res_.resize(domain_.mesh().ncells(), 0.0);
        
        // Create fields with ghost cells
        u_ = domain_.create_field_with_ghosts<Vector>();
        un_ = domain_.create_field_with_ghosts<Vector>();
    }
    
    // ===================================================================
    // Getters
    // ===================================================================
    const ComputationalDomain& domain() const { return domain_; }
    
    Vector& q_face_left() { return q_face_left_; }
    Vector& q_face_right() { return q_face_right_; }
    Vector& flux() { return flux_; }
    Vector& res() { return res_; }
    Vector& u() { return u_; }
    Vector& un() { return un_; }
    
    const Vector& q_face_left() const { return q_face_left_; }
    const Vector& q_face_right() const { return q_face_right_; }
    const Vector& flux() const { return flux_; }
    const Vector& res() const { return res_; }
    const Vector& u() const { return u_; }
    const Vector& un() const { return un_; }
    
    // ===================================================================
    // Solution management
    // ===================================================================
    
    // Copy current solution to old solution
    void copy_to_old() {
        std::copy(u_.begin(), u_.end(), un_.begin());
    }
    
    // Swap current and old solutions
    void swap_solutions() {
        std::swap(u_, un_);
    }
    
    // Get interior solution (without ghost cells)
    Vector get_interior_solution() const {
        return domain_.extract_interior(u_);
    }
    
    // Get cell center coordinates
    const Vector& get_cell_centers() const {
        return domain_.mesh().cell_centers();
    }
    
    // ===================================================================
    // Boundary conditions
    // ===================================================================
    
    // Apply periodic boundary conditions
    void apply_periodic_boundary() {
        domain_.apply_periodic_boundary(u_);
    }
    
    // Apply boundary conditions (wrapper)
    void apply_boundary() {
        apply_periodic_boundary();
    }
    
    // ===================================================================
    // Initial conditions
    // ===================================================================
    
    void initialize(const std::unique_ptr<InitialCondition>& ic,
                   const std::unique_ptr<BoundaryCondition>& bc);

    // ===================================================================
    // Validation and diagnostics
    // ===================================================================
    
    // Validate solution
    void validate() const {
        if (u_.size() != static_cast<size_t>(domain_.ntcells())) {
            throw std::runtime_error("Solution size does not match domain");
        }
        
        if (q_face_left_.size() != static_cast<size_t>(domain_.mesh().nnodes())) {
            throw std::runtime_error("Face left array size mismatch");
        }
        
        if (q_face_right_.size() != static_cast<size_t>(domain_.mesh().nnodes())) {
            throw std::runtime_error("Face right array size mismatch");
        }
        
        if (flux_.size() != static_cast<size_t>(domain_.mesh().nnodes())) {
            throw std::runtime_error("Flux array size mismatch");
        }
        
        if (res_.size() != static_cast<size_t>(domain_.mesh().ncells())) {
            throw std::runtime_error("Residual array size mismatch");
        }
    }
    
    // Print solution information
    void print_info() const {
        std::cout << "Solution Information:" << std::endl;
        std::cout << "  Field size (with ghosts): " << u_.size() << std::endl;
        std::cout << "  Interior cells: " << domain_.mesh().ncells() << std::endl;
        std::cout << "  Ghost cells: " << domain_.nghosts() << " each side" << std::endl;
        
        // Compute statistics
        auto interior = get_interior_solution();
        auto [min_val, max_val] = std::minmax_element(interior.begin(), interior.end());
        Real mean = std::accumulate(interior.begin(), interior.end(), 0.0) / interior.size();
        
        std::cout << "  Interior solution statistics:" << std::endl;
        std::cout << "    Min: " << *min_val << std::endl;
        std::cout << "    Max: " << *max_val << std::endl;
        std::cout << "    Mean: " << mean << std::endl;
    }
    
    // Print detailed solution values
    void print_detailed(int max_cells = 10) const {
        const auto& mesh = domain_.mesh();
        auto interior = get_interior_solution();
        
        std::cout << "Detailed Solution:" << std::endl;
        std::cout << "  Index | X-coordinate | Solution Value" << std::endl;
        std::cout << "  -------------------------------------" << std::endl;
        
        int n_to_print = std::min(max_cells, mesh.ncells());
        for (int i = 0; i < n_to_print; ++i) {
            std::cout << "  " << std::setw(5) << i 
                      << " | " << std::setw(12) << std::setprecision(6) << mesh.cell_center_coordinate(i)
                      << " | " << std::setw(12) << std::setprecision(6) << interior[i] 
                      << std::endl;
        }
        
        if (mesh.ncells() > n_to_print) {
            std::cout << "  ... (showing first " << n_to_print << " of " << mesh.ncells() << " cells)" << std::endl;
        }
    }
    
    // ===================================================================
    // File I/O
    // ===================================================================
    
    // Write solution to file
    void write_to_file(const std::string& filename, 
                      const Vector& solution,
                      const std::string& comment = "") const {
        std::ofstream file(filename);
        
        if (!comment.empty()) {
            file << "# " << comment << std::endl;
        }
        
        file << "# Solution Data" << std::endl;
        file << "# Domain: [" << domain_.mesh().xmin() << ", " << domain_.mesh().xmax() << "]" << std::endl;
        file << "# Cells: " << domain_.mesh().ncells() << std::endl;
        file << "# dx: " << domain_.mesh().dx() << std::endl;
        file << "#" << std::endl;
        //file << "# Index, X-coordinate, Solution" << std::endl;
        file << "# X-coordinate, Solution" << std::endl;
        
        file << std::fixed << std::setprecision(6);
        const auto& xcc = domain_.mesh().cell_centers();
        
        for (int i = 0; i < domain_.mesh().ncells(); ++i) {
            //file << std::setw(6) << i 
            file << std::setw(12) << xcc[i]
                 << std::setw(12) << solution[i] 
                 << std::endl;
        }
        
        file.close();
    }
    
    // Write current interior solution to file
    void write_current_solution(const std::string& filename,
                               const std::string& comment = "") const {
        auto interior = get_interior_solution();
        write_to_file(filename, interior, comment);
    }
    
    // ===================================================================
    // Factory methods
    // ===================================================================
    static Solution create_from_domain(const ComputationalDomain& domain) {
        return Solution(domain);
    }
};

} // namespace cfd

#endif // SOLUTION_HPP