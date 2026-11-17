#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include "cfd_solver.hpp"
#include <vector>
#include <string>
#include <memory>

namespace cfd {

// ===================================================================
// Analysis utilities for CFD simulations
// ===================================================================
class Analysis {
public:
    // ===================================================================
    // Run ENO vs WENO comparison analysis
    // ===================================================================
    static void run_eno_weno(int cells = 40, double final_time = 0.625, 
                            bool verbose = true);
							
	static void run_complex_wave_analysis(int cells, double final_time, bool verbose);
    
    // ===================================================================
    // Run convergence study
    // ===================================================================
    static void run_convergence_study(const std::string& scheme,
                                     const std::vector<int>& cell_counts = {20, 40, 80, 160},
                                     double final_time = 0.625);
    
    // ===================================================================
    // Run flux comparison
    // ===================================================================
    static void run_flux_comparison(const std::string& recon_scheme,
                                   int cells = 40,
                                   double final_time = 0.625);
    
    // ===================================================================
    // Run time integrator comparison
    // ===================================================================
    static void run_time_integrator_comparison(const std::string& recon_scheme,
                                              int cells = 40,
                                              double final_time = 0.625);
    
    // ===================================================================
    // Run sensitivity analysis
    // ===================================================================
    static void run_sensitivity_analysis(const CfdConfig& base_config,
                                        const std::string& parameter,
                                        const std::vector<double>& values);
    
    // ===================================================================
    // Utility functions
    // ===================================================================
    static double compute_error(const Vector& numerical, const Vector& analytical);
    static double compute_l1_error(const Vector& numerical, const Vector& analytical);
    static double compute_l2_error(const Vector& numerical, const Vector& analytical);
    static double compute_linf_error(const Vector& numerical, const Vector& analytical);
    
    // ===================================================================
    // File I/O for analysis results
    // ===================================================================
    static void write_analysis_results(const std::string& filename,
                                      const std::vector<double>& errors,
                                      const std::vector<int>& cells);
    
    static void plot_results_script(const std::string& output_dir = "results");
};

} // namespace cfd

#endif // ANALYSIS_HPP