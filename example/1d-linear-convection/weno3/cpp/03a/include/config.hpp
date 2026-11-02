// ==================== include/config.hpp ====================
#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <algorithm>

namespace cfd {

// Type aliases
using Real = double;
using Vector = std::vector<Real>;

// ===================================================================
// Configuration structure for CFD solver
// ===================================================================
struct CfdConfig {
    // Initial condition type
    std::string ic_type = "complex";  // "step", "sine", "complex", etc.
    
    // Reconstruction scheme
    std::string recon_scheme = "eno";
    
    // Flux type - 改为字符串
    std::string flux_type = "rusanov";  // "rusanov", "engquist-osher", "upwind", "central", "roe", "hll"
    
    // Time integration method - 新增
    std::string time_integration = "rk2";  // "rk1", "rk2", "rk3", "crank-nicolson", "implicit-euler"
    
    // Theta method parameter (0.0=explicit, 0.5=C-N, 1.0=implicit)
    Real theta = 0.5;
    
    // Time integration order (保持向后兼容)
    int rk_order = 2;
    
    // Spatial order of accuracy
    int spatial_order = 3;
    
    // Wave speed (advection velocity)
    double wave_speed = 1.0;
    
    // Simulation final time
    double final_time = 0.0;
    
    // Time step size
    double dt = 0.025;
    
    // CFL number for stability
    double cfl = 0.5;
    
    // Mesh parameters
    double xmin = 0.0;
    double xmax = 2.0;
    int ncells = 40;
    
    // Output control
    bool verbose = false;
    int print_interval = 100;
    
    // File output options
    std::string output_dir = ".";
    std::string eno_filename = "eno_results.txt";
    std::string weno_filename = "weno_results.txt";
    std::string analytical_filename = "analytical_results.txt";
    std::string plot_filename = "comparison.png";
    
    // 初始条件参数（可选）
    std::vector<Real> ic_params;
    
    // ===================================================================
    // Constructor
    // ===================================================================
    CfdConfig() = default;
    
    // Constructor with basic parameters
    CfdConfig(const std::string& ic_type, 
              const std::string& scheme, 
              int cells, 
              double time) 
        : ic_type(ic_type), recon_scheme(scheme), ncells(cells), final_time(time) {}
    
    // ===================================================================
    // Validation methods
    // ===================================================================
    void validate() const {
        // Validate initial condition type
        std::vector<std::string> valid_ic_types = {
            "step", "sine", "sinusoidal", "gaussian", 
            "complex", "complexwave", "custom"
        };
        
        bool valid_ic = false;
        std::string ic_lower = ic_type;
        std::transform(ic_lower.begin(), ic_lower.end(), ic_lower.begin(), ::tolower);
        
        for (const auto& valid_type : valid_ic_types) {
            if (ic_lower == valid_type) {
                valid_ic = true;
                break;
            }
        }
        
        if (!valid_ic) {
            throw std::runtime_error("Invalid initial condition type: " + ic_type + 
                                   "\nValid types: step, sine/sinusoidal, gaussian, complex/complexwave, custom");
        }
        
        // Validate reconstruction scheme
        if (recon_scheme != "eno" && recon_scheme != "weno") {
            throw std::runtime_error("Invalid reconstruction scheme: " + recon_scheme + 
                                   " (must be 'eno' or 'weno')");
        }
        
        // Validate flux type (新验证)
        std::vector<std::string> valid_flux_types = {
            "rusanov", "engquist", "engquist-osher", "eo",
            "upwind", "central", "lax-wendroff", "lw",
            "roe", "hll"
        };
        
        std::string flux_lower = flux_type;
        std::transform(flux_lower.begin(), flux_lower.end(), flux_lower.begin(), ::tolower);
        
        bool valid_flux = false;
        for (const auto& valid_type : valid_flux_types) {
            if (flux_lower == valid_type) {
                valid_flux = true;
                break;
            }
        }
        
        if (!valid_flux) {
            throw std::runtime_error("Invalid flux type: " + flux_type + 
                                   "\nValid types: rusanov, engquist-osher, upwind, central, lax-wendroff, roe, hll");
        }
        
        // Validate time integration method
        std::vector<std::string> valid_time_methods = {
            "rk1", "rk2", "rk3", "crank-nicolson", "cn",
            "implicit-euler", "explicit-euler", "theta",
            "euler", "forward-euler", "midpoint", "tvdrk3"
        };
        
        std::string time_lower = time_integration;
        std::transform(time_lower.begin(), time_lower.end(), time_lower.begin(), ::tolower);
        
        bool valid_time = false;
        for (const auto& valid_method : valid_time_methods) {
            if (time_lower == valid_method) {
                valid_time = true;
                break;
            }
        }
        
        if (!valid_time) {
            throw std::runtime_error("Invalid time integration method: " + time_integration + 
                                   "\nValid methods: rk1, rk2, rk3, crank-nicolson, implicit-euler, explicit-euler, theta");
        }
        
        // Validate spatial order
        if (recon_scheme == "eno") {
            // ENO 支持 1, 2, 3, 5 阶
            if (spatial_order != 1 && spatial_order != 2 && spatial_order != 3 && spatial_order != 5) {
                throw std::runtime_error("ENO spatial order must be 1, 2, 3, or 5, got: " +
                    std::to_string(spatial_order));
            }
        } else if (recon_scheme == "weno") {
            // WENO 支持 1, 3, 5 阶
            if (spatial_order != 1 && spatial_order != 3 && spatial_order != 5) {
                throw std::runtime_error("WENO spatial order must be 1, 3, or 5, got: " +
                    std::to_string(spatial_order));
            }
        }
        
        // Validate mesh parameters
        if (xmin >= xmax) {
            throw std::runtime_error("xmin must be less than xmax");
        }
        
        if (ncells <= 0) {
            throw std::runtime_error("Number of cells must be positive, got: " +
                                   std::to_string(ncells));
        }
        
        // Validate time parameters
        if (final_time < 0.0) {
            throw std::runtime_error("Final time cannot be negative, got: " +
                                   std::to_string(final_time));
        }
        
        if (dt <= 0.0) {
            throw std::runtime_error("Time step must be positive, got: " +
                                   std::to_string(dt));
        }
        
        if (cfl <= 0.0 || cfl > 1.0) {
            throw std::runtime_error("CFL number must be in (0, 1], got: " +
                                   std::to_string(cfl));
        }
        
        // Validate wave speed
        if (wave_speed == 0.0) {
            throw std::runtime_error("Wave speed cannot be zero");
        }
        
        // Validate theta parameter
        if (theta < 0.0 || theta > 1.0) {
            throw std::runtime_error("Theta parameter must be in [0, 1], got: " +
                                   std::to_string(theta));
        }
    }
    
    // ===================================================================
    // Helper methods for flux type conversion (保持向后兼容)
    // ===================================================================
    int flux_type_to_int() const {
        std::string flux_lower = flux_type;
        std::transform(flux_lower.begin(), flux_lower.end(), flux_lower.begin(), ::tolower);
        
        if (flux_lower == "rusanov" || flux_lower == "llf") {
            return 0;
        } else if (flux_lower == "engquist" || flux_lower == "engquist-osher" || flux_lower == "eo") {
            return 1;
        } else if (flux_lower == "lax-wendroff" || flux_lower == "lw") {
            return 2;
        } else if (flux_lower == "upwind") {
            return 3;
        } else if (flux_lower == "central") {
            return 4;
        } else if (flux_lower == "roe") {
            return 5;
        } else if (flux_lower == "hll") {
            return 6;
        } else {
            return 0;  // Default to Rusanov
        }
    }
    
    static std::string flux_int_to_string(int flux_int) {
        switch (flux_int) {
            case 0: return "rusanov";
            case 1: return "engquist-osher";
            case 2: return "lax-wendroff";
            case 3: return "upwind";
            case 4: return "central";
            case 5: return "roe";
            case 6: return "hll";
            default: return "rusanov";
        }
    }
    
    // ===================================================================
    // Helper methods for common configurations
    // ===================================================================
    
    // Create step function configuration
    static CfdConfig create_step_config(int cells = 40, double final_time = 0.625) {
        CfdConfig config;
        config.ic_type = "step";
        config.recon_scheme = "eno";
        config.spatial_order = 3;
        config.flux_type = "rusanov";
        config.time_integration = "rk2";
        config.rk_order = 2;
        config.wave_speed = 1.0;
        config.final_time = final_time;
        config.cfl = 1.0;
        config.dt = 0.0025;
        config.ncells = cells;
        config.xmin = 0.0;
        config.xmax = 2.0;
        
        // Step function parameters: low_val, high_val, start, end
        config.ic_params = {1.0, 2.0, 0.5, 1.0};
        
        return config;
    }
    
    // Create Crank-Nicolson configuration
    static CfdConfig create_crank_nicolson_config(int cells = 40, double final_time = 0.625) {
        CfdConfig config = create_step_config(cells, final_time);
        config.time_integration = "crank-nicolson";
        config.flux_type = "central";  // C-N 通常与中心差分配合使用
        config.dt = 0.01;  // 可以使用更大的时间步长
        return config;
    }
    
    // Create ENO configuration
    static CfdConfig create_eno_config(int order = 3, int cells = 40, 
                                      double final_time = 0.625) {
        CfdConfig config = create_step_config(cells, final_time);
        config.recon_scheme = "eno";
        config.spatial_order = order;
        config.rk_order = 2;
        config.time_integration = "rk2";
        return config;
    }
    
    // Create WENO configuration
    static CfdConfig create_weno_config(int cells = 40, double final_time = 0.625) {
        CfdConfig config = create_step_config(cells, final_time);
        config.recon_scheme = "weno";
        config.spatial_order = 3;
        config.rk_order = 2;
        config.time_integration = "rk2";
        return config;
    }
    
    // ===================================================================
    // Print methods
    // ===================================================================
    void print() const {
        std::cout << "==========================================" << std::endl;
        std::cout << "CFD Solver Configuration" << std::endl;
        std::cout << "==========================================" << std::endl;
        std::cout << "  Initial Condition: " << ic_type << std::endl;
        std::cout << "  Reconstruction: " << recon_scheme << std::endl;
        std::cout << "  Spatial Order: " << spatial_order << std::endl;
        std::cout << "  Flux Type: " << flux_type << std::endl;
        std::cout << "  Time Integration: " << time_integration << std::endl;
        if (time_integration.find("theta") != std::string::npos || 
            time_integration == "crank-nicolson") {
            std::cout << "  Theta Parameter: " << theta << std::endl;
        }
        std::cout << "  Wave Speed: " << wave_speed << std::endl;
        std::cout << "  Final Time: " << final_time << std::endl;
        std::cout << "  Time Step: " << dt << std::endl;
        std::cout << "  CFL Number: " << cfl << std::endl;
        std::cout << "  Mesh Cells: " << ncells << std::endl;
        std::cout << "  Domain: [" << xmin << ", " << xmax << "]" << std::endl;
        std::cout << "  Domain Length: " << domain_length() << std::endl;
        
        if (!ic_params.empty()) {
            std::cout << "  IC Parameters: [";
            for (size_t i = 0; i < ic_params.size(); ++i) {
                std::cout << ic_params[i];
                if (i < ic_params.size() - 1) std::cout << ", ";
            }
            std::cout << "]" << std::endl;
        }
        
        if (verbose) {
            std::cout << "  Verbose Output: Yes" << std::endl;
            std::cout << "  Print Interval: " << print_interval << " steps" << std::endl;
        }
        std::cout << "==========================================" << std::endl;
    }

    void print_brief() const {
        std::cout << ic_type << " + " << recon_scheme << spatial_order 
            << " (cells=" << ncells 
            << ", t=" << final_time
            << ", dt=" << dt 
            << ", CFL=" << cfl << ")";
    }
    
    // ===================================================================
    // Getters for derived parameters
    // ===================================================================
    double domain_length() const { return xmax - xmin; }
    double cell_size() const { return domain_length() / ncells; }
    
    // Check if configuration is valid
    bool is_valid() const {
        try {
            validate();
            return true;
        } catch (const std::exception&) {
            return false;
        }
    }
    
    // Check if time integration is implicit
    bool is_implicit() const {
        std::string time_lower = time_integration;
        std::transform(time_lower.begin(), time_lower.end(), time_lower.begin(), ::tolower);
        
        return (time_lower == "crank-nicolson" || 
                time_lower == "cn" ||
                time_lower == "implicit-euler" ||
                time_lower == "theta" ||
                time_lower.find("theta-") == 0);
    }
    
    // Check if this is a complex wave configuration
    bool is_complex_wave() const {
        std::string ic_lower = ic_type;
        std::transform(ic_lower.begin(), ic_lower.end(), ic_lower.begin(), ::tolower);
        return (ic_lower == "complex" || ic_lower == "complexwave");
    }
};

} // namespace cfd

#endif // CONFIG_HPP