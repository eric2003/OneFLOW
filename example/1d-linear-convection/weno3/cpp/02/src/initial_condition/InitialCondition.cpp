// ==================== src/initial_condition/InitialCondition.cpp ====================
#include "initial_condition/InitialCondition.h"

namespace cfd {

// ===================================================================
// Default implementation for initialize_with_ghosts
// ===================================================================
void InitialCondition::initialize_with_ghosts(Vector& u_with_ghosts, 
                                             const ComputationalDomain& domain) const {
    // Check if the vector has the correct size
    if (u_with_ghosts.size() != static_cast<size_t>(domain.ntcells())) {
        u_with_ghosts.resize(domain.ntcells(), 0.0);
    }
    
    // Get interior portion of the array
    const int ist = domain.ist();
    const int ied = domain.ied();
    
    // Extract interior indices
    Vector u_interior(domain.interior_cells());
    
    // Initialize interior solution using the derived class implementation
    initialize(u_interior, domain);
    
    // Copy interior values to the ghosted array
    for (int i = 0; i < domain.interior_cells(); ++i) {
        u_with_ghosts[ist + i] = u_interior[i];
    }
    
    // Apply boundary conditions (ghost cells will be filled by boundary conditions later)
    // Note: The ghost cells remain uninitialized here; they should be filled
    // by the boundary condition class during the simulation
}

} // namespace cfd