#ifndef GAUSSIAN_INITIAL_CONDITION_H
#define GAUSSIAN_INITIAL_CONDITION_H

#include "InitialCondition.h"

namespace cfd {

// ===================================================================
// Gaussian initial condition
// ===================================================================
class GaussianInitialCondition : public InitialCondition {
private:
    Real amplitude_;
    Real center_;
    Real width_;
    Real background_;
    
public:
    GaussianInitialCondition(Real amplitude = 1.0, Real center = 1.0,
                            Real width = 0.1, Real background = 1.0);
    
    ~GaussianInitialCondition() override = default;
    
    void initialize(Vector& u_interior, const ComputationalDomain& domain) const override;
    
    Real evaluate(Real x) const override;
    
    std::string name() const override;
    
    int type_id() const override { return 2; }
    
    std::unique_ptr<InitialCondition> clone() const override;
    
    // Getters
    Real amplitude() const { return amplitude_; }
    Real center() const { return center_; }
    Real width() const { return width_; }
    Real background() const { return background_; }
    
    // Setters
    void set_amplitude(Real val) { amplitude_ = val; }
    void set_center(Real val) { center_ = val; }
    void set_width(Real val) { width_ = val; }
    void set_background(Real val) { background_ = val; }
};

} // namespace cfd

#endif // GAUSSIAN_INITIAL_CONDITION_H