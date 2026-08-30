// ==================== src/solution.cpp ====================
#include "solution.hpp"
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <fstream>

namespace cfd {

void Solution::initialize(const std::unique_ptr<InitialCondition>& ic,
                         const std::unique_ptr<BoundaryCondition>& bc) {
    if (ic) {
        ic->initialize_with_ghosts(u_, domain_);
    }
    
    if (bc) {
        bc->apply(u_, domain_);
    }
    
    copy_to_old();
}

} // namespace cfd