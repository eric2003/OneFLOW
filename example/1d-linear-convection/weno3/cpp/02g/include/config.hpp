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
    
    // Flux type: 0 = Rusanov, 1 = Engquist-Osher
    int flux_type = 0;
    
    // Time integration order
    int rk_order = 1;
    
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
    
    // ===================================================================
    // 初始条件参数（可选）
    // ===================================================================
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
        
        // Validate flux type
        if (flux_type != 0 && flux_type != 1) {
            throw std::runtime_error("Invalid flux type: " + std::to_string(flux_type) +
                                   " (must be 0 or 1)");
        }
        
        // Validate RK order
        if (rk_order != 1 && rk_order != 2 && rk_order != 3 ) {
            throw std::runtime_error("Invalid RK order: " + std::to_string(rk_order) +
                                   " (must be 1 or 2)");
        }
        
        // Validate spatial order
        if (recon_scheme == "eno") {
            if (spatial_order < 1 || spatial_order > 3) {
                throw std::runtime_error("ENO spatial order must be 1, 2, or 3, got: " +
                                       std::to_string(spatial_order));
            }
        } else if (recon_scheme == "weno") {
            if (spatial_order != 3) {
                throw std::runtime_error("WENO spatial order must be 3, got: " +
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
        config.flux_type = 0;
        config.rk_order = 1;
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
    
    // Create sine wave configuration
    static CfdConfig create_sine_config(int cells = 40, double final_time = 0.625) {
        CfdConfig config;
        config.ic_type = "sine";
        config.recon_scheme = "eno";
        config.spatial_order = 3;
        config.flux_type = 0;
        config.rk_order = 1;
        config.wave_speed = 1.0;
        config.final_time = final_time;
        config.cfl = 1.0;
        config.dt = 0.0025;
        config.ncells = cells;
        config.xmin = 0.0;
        config.xmax = 2.0;
        
        // Sine wave parameters: amplitude, frequency, mean, phase
        config.ic_params = {0.5, 1.0, 1.0, 0.0};
        
        return config;
    }
    
    // Create complex wave configuration (from Fortran)
    static CfdConfig create_complex_config(int cells = 40, double final_time = 0.0) {
        CfdConfig config;
        config.ic_type = "complex";
        config.recon_scheme = "eno";
        config.spatial_order = 3;
        config.flux_type = 0;
        config.rk_order = 1;
        config.wave_speed = 1.0;
        config.final_time = final_time;
        config.cfl = 1.0;
        config.dt = 0.0025;
        config.ncells = cells;
        config.xmin = -1.0;  // Complex wave needs [-1, 1] domain
        config.xmax = 1.0;
        
        // Complex wave parameters: a, z, delta, alpha
        // a = 0.5, z = -0.7, delta = 0.005, alpha = 10.0
        config.ic_params = {0.5, -0.7, 0.005, 10.0};
        
        return config;
    }
    
    // Create ENO configuration
    static CfdConfig create_eno_config(int order = 3, int cells = 40, 
                                      double final_time = 0.625) {
        CfdConfig config = create_step_config(cells, final_time);
        config.recon_scheme = "eno";
        config.spatial_order = order;
        config.rk_order = 2;
        return config;
    }
    
    // Create WENO configuration
    static CfdConfig create_weno_config(int cells = 40, double final_time = 0.625) {
        CfdConfig config = create_step_config(cells, final_time);
        config.recon_scheme = "weno";
        config.spatial_order = 3;
        config.rk_order = 2;
        return config;
    }
    
    // Create configuration for analysis (ENO vs WENO comparison)
    static std::pair<CfdConfig, CfdConfig> create_analysis_configs(
        int cells = 40, double final_time = 0.625) {
        return std::make_pair(
            create_eno_config(3, cells, final_time),
            create_weno_config(cells, final_time)
        );
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
        std::cout << "  Flux Type: " << (flux_type == 0 ? "Rusanov" : "Engquist-Osher") << std::endl;
        std::cout << "  RK Order: " << rk_order << std::endl;
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
    
    // Check if this is a complex wave configuration
    bool is_complex_wave() const {
        std::string ic_lower = ic_type;
        std::transform(ic_lower.begin(), ic_lower.end(), ic_lower.begin(), ::tolower);
        return (ic_lower == "complex" || ic_lower == "complexwave");
    }
};

} // namespace cfd

#endif // CONFIG_HPP