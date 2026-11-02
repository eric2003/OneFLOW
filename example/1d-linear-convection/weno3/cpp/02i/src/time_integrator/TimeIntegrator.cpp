#include "time_integrator/TimeIntegrator.h"

namespace cfd {

TimeIntegrator::TimeIntegrator(std::unique_ptr<ResidualCalculator> residual_calculator)
    : residual_calculator_(std::move(residual_calculator)) {}

const ResidualCalculator& TimeIntegrator::residual_calculator() const { 
    return *residual_calculator_; 
}

} // namespace cfd