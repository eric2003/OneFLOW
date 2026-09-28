#ifndef RECONSTRUCTOR_FACTORY_H
#define RECONSTRUCTOR_FACTORY_H

#include "Reconstructor.h"
#include "basic/FirstOrderReconstructor.h"
#include "basic/SecondOrderReconstructor.h"
#include "basic/MUSCLReconstructor.h"
#include "advanced/EnoReconstructor.h"
#include "advanced/WenoReconstructor.h"
#include "advanced/MPWenoReconstructor.h"
#include <memory>
#include <vector>
#include <string>

namespace cfd {

// ===================================================================
// Reconstructor factory
// ===================================================================
class ReconstructorFactory {
public:
    // Create reconstructor by name with parameters
    static std::unique_ptr<Reconstructor> create_reconstructor(
        const std::string& scheme, 
        int order = 3, 
        const std::string& limiter = "minmod",
        Real param = 1.0);
    
    // Create reconstructor from configuration
    static std::unique_ptr<Reconstructor> create_reconstructor(
        const CfdConfig& config);
    
    // Create reconstructor by order (simple interface)
    static std::unique_ptr<Reconstructor> create_by_order(
        int order, 
        const std::string& type = "linear");
    
    // Get available reconstructors
    static std::vector<std::string> available_reconstructors();
    static std::vector<std::string> available_reconstructors_by_type(const std::string& type);
    
    // Get reconstructor information
    static std::string get_reconstructor_info(const std::string& scheme, int order = 3);
    
    // Check if reconstructor is available
    static bool is_available(const std::string& scheme, int order = -1);
    
    // Get default reconstructor
    static std::unique_ptr<Reconstructor> create_default_reconstructor();
    
    // Create reconstructor for specific purpose
    static std::unique_ptr<Reconstructor> create_for_shocks();
    static std::unique_ptr<Reconstructor> create_for_smooth();
    static std::unique_ptr<Reconstructor> create_high_order();
    static std::unique_ptr<Reconstructor> create_tvd();
    
    // Comparison utilities
    static void compare_reconstructors(const Vector& q,
                                      const ComputationalDomain& domain,
                                      const std::vector<std::string>& schemes);
};

} // namespace cfd

#endif // RECONSTRUCTOR_FACTORY_H