#include "reconstructor/ReconstructorFactory.h"
#include <algorithm>
#include <iostream>
#include <sstream>

namespace cfd {

std::unique_ptr<Reconstructor> ReconstructorFactory::create_reconstructor(
    const std::string& scheme, 
    int order, 
    const std::string& limiter,
    Real param) {
    
    std::string scheme_lower = scheme;
    std::transform(scheme_lower.begin(), scheme_lower.end(), 
                   scheme_lower.begin(), ::tolower);
    
    if (scheme_lower == "first" || scheme_lower == "first-order" || scheme_lower == "1") {
        return std::make_unique<FirstOrderReconstructor>();
    }
    else if (scheme_lower == "second" || scheme_lower == "second-order" || scheme_lower == "2") {
        return std::make_unique<SecondOrderReconstructor>(param, limiter != "none");
    }
    else if (scheme_lower == "muscl") {
        return std::make_unique<MUSCLReconstructor>(order, limiter);
    }
    else if (scheme_lower == "eno") {
        return std::make_unique<EnoReconstructor>(order);
    }
    else if (scheme_lower == "weno") {
        return std::make_unique<WenoReconstructor>(order);
    }
    else if (scheme_lower == "mpweno") {
        return std::make_unique<MPWenoReconstructor>(order);
    }
    else {
        throw std::invalid_argument("Unknown reconstruction scheme: " + scheme);
    }
}

std::unique_ptr<Reconstructor> ReconstructorFactory::create_reconstructor(
    const CfdConfig& config) {
    
    return create_reconstructor(
        config.recon_scheme,
        config.spatial_order,
        "minmod",  // default limiter
        1.0        // default parameter
    );
}

std::unique_ptr<Reconstructor> ReconstructorFactory::create_by_order(
    int order, 
    const std::string& type) {
    
    if (order == 1) {
        return std::make_unique<FirstOrderReconstructor>();
    } 
    else if (order == 2) {
        if (type == "linear") {
            return std::make_unique<SecondOrderReconstructor>();
        } else if (type == "tvd") {
            return std::make_unique<MUSCLReconstructor>(2, "minmod");
        }
    }
    else if (order == 3) {
        if (type == "linear") {
            return std::make_unique<EnoReconstructor>(3);
        } else if (type == "nonlinear") {
            return std::make_unique<WenoReconstructor>(3);
        } else if (type == "mp") {
            return std::make_unique<MPWenoReconstructor>(3);
        }
    }
    
    throw std::invalid_argument("Unsupported order or type: order=" + 
                               std::to_string(order) + ", type=" + type);
}

std::vector<std::string> ReconstructorFactory::available_reconstructors() {
    return {
        "First-Order (Piecewise Constant)",
        "Second-Order (Linear)",
        "MUSCL (TVD)",
        "ENO3 (Essentially Non-Oscillatory)",
        "WENO3 (Weighted ENO)",
        "MP-WENO3 (Monotonicity Preserving WENO)"
    };
}

std::vector<std::string> ReconstructorFactory::available_reconstructors_by_type(const std::string& type) {
    std::vector<std::string> all = {
        "first", "second", "muscl", "eno", "weno", "mpweno"
    };
    
    std::vector<std::string> result;
    
    for (const auto& scheme : all) {
        try {
            auto recon = create_reconstructor(scheme, 3, "minmod");
            if (type == "linear" && recon->is_linear()) {
                result.push_back(scheme);
            } else if (type == "nonlinear" && recon->is_nonlinear()) {
                result.push_back(scheme);
            } else if (type == "tvd" && recon->is_tvd()) {
                result.push_back(scheme);
            } else if (type == "high-order" && recon->order() >= 3) {
                result.push_back(scheme);
            }
        } catch (...) {
            continue;
        }
    }
    
    return result;
}

std::string ReconstructorFactory::get_reconstructor_info(const std::string& scheme, int order) {
    try {
        auto recon = create_reconstructor(scheme, order);
        std::ostringstream oss;
        oss << recon->name() << "\n";
        oss << "  Order: " << recon->order() << "\n";
        oss << "  Stencil width: " << recon->stencil_width() << "\n";
        oss << "  Ghost cells required: " << recon->required_ghost_cells() << "\n";
        oss << "  Linear: " << (recon->is_linear() ? "Yes" : "No") << "\n";
        oss << "  TVD: " << (recon->is_tvd() ? "Yes" : "No") << "\n";
        oss << "  Monotonic: " << (recon->is_monotonic() ? "Yes" : "No");
        return oss.str();
    } catch (const std::exception& e) {
        return "Error: " + std::string(e.what());
    }
}

bool ReconstructorFactory::is_available(const std::string& scheme, int order) {
    try {
        auto recon = create_reconstructor(scheme, order);
        return recon != nullptr;
    } catch (...) {
        return false;
    }
}

std::unique_ptr<Reconstructor> ReconstructorFactory::create_default_reconstructor() {
    return std::make_unique<WenoReconstructor>(3);
}

std::unique_ptr<Reconstructor> ReconstructorFactory::create_for_shocks() {
    return std::make_unique<MPWenoReconstructor>(3);
}

std::unique_ptr<Reconstructor> ReconstructorFactory::create_for_smooth() {
    return std::make_unique<EnoReconstructor>(3);
}

std::unique_ptr<Reconstructor> ReconstructorFactory::create_high_order() {
    return std::make_unique<WenoReconstructor>(3);
}

std::unique_ptr<Reconstructor> ReconstructorFactory::create_tvd() {
    return std::make_unique<MUSCLReconstructor>(2, "minmod");
}

void ReconstructorFactory::compare_reconstructors(const Vector& q,
                                                const ComputationalDomain& domain,
                                                const std::vector<std::string>& schemes) {
    std::cout << "\nReconstructor Comparison" << std::endl;
    std::cout << "========================" << std::endl;
    
    for (const auto& scheme : schemes) {
        try {
            auto recon = create_reconstructor(scheme, 3);
            Vector q_left(q.size()), q_right(q.size());
            recon->reconstruct(q, q_left, q_right, domain);
            
            std::cout << scheme << ": ";
            recon->print_info();
        } catch (const std::exception& e) {
            std::cout << scheme << ": ERROR - " << e.what() << std::endl;
        }
    }
}

} // namespace cfd