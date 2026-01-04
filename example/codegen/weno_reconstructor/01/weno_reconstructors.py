# Auto-generated WENO coefficient initializers
# DO NOT EDIT MANUALLY

def _reconstruct_weno1_left(self, v1):
    eps = 1e-6
    beta0 = 0
    d0 = 1.0
    alpha0 = d0 / (eps + beta0)**2
    alpha = alpha0
    w0 = alpha0 / alpha
    q0 = 1.0*v1  # r=0
    return w0 * q0

def _reconstruct_weno1_right(self, v1):
    eps = 1e-6
    beta0 = 0
    d0 = 1.0
    alpha0 = d0 / (eps + beta0)**2
    alpha = alpha0
    w0 = alpha0 / alpha
    q0 = 1.0*v1  # r=0
    return w0 * q0

def _reconstruct_weno3_left(self, v1, v2, v3):
    eps = 1e-6
    beta0 = (v2 - v1)**2
    beta1 = (v3 - v2)**2
    d0 = 1.0/3.0
    d1 = 2.0/3.0
    alpha0 = d0 / (eps + beta0)**2
    alpha1 = d1 / (eps + beta1)**2
    alpha = alpha0 + alpha1
    w0 = alpha0 / alpha
    w1 = alpha1 / alpha
    q0 = -1.0/2.0*v1+3.0/2.0*v2  # r=1
    q1 = 1.0/2.0*v2+1.0/2.0*v3  # r=0
    return w0 * q0 + w1 * q1

def _reconstruct_weno3_right(self, v1, v2, v3):
    eps = 1e-6
    beta0 = (v2 - v1)**2
    beta1 = (v3 - v2)**2
    d0 = 2.0/3.0
    d1 = 1.0/3.0
    alpha0 = d0 / (eps + beta0)**2
    alpha1 = d1 / (eps + beta1)**2
    alpha = alpha0 + alpha1
    w0 = alpha0 / alpha
    w1 = alpha1 / alpha
    q0 = 1.0/2.0*v1+1.0/2.0*v2  # r=1
    q1 = 3.0/2.0*v2-1.0/2.0*v3  # r=0
    return w0 * q0 + w1 * q1

def _reconstruct_weno5_left(self, v1, v2, v3, v4, v5):
    eps = 1e-6
    beta0 = (13.0/12.0)*(v1 - 2*v2 + v3)**2 + (1.0/4.0)*(v1 - 4*v2 + 3*v3)**2
    beta1 = (13.0/12.0)*(v2 - 2*v3 + v4)**2 + (1.0/4.0)*(v2 - v4)**2
    beta2 = (13.0/12.0)*(v3 - 2*v4 + v5)**2 + (1.0/4.0)*(3*v3 - 4*v4 + v5)**2
    d0 = 1.0/10.0
    d1 = 3.0/5.0
    d2 = 3.0/10.0
    alpha0 = d0 / (eps + beta0)**2
    alpha1 = d1 / (eps + beta1)**2
    alpha2 = d2 / (eps + beta2)**2
    alpha = alpha0 + alpha1 + alpha2
    w0 = alpha0 / alpha
    w1 = alpha1 / alpha
    w2 = alpha2 / alpha
    q0 = 1.0/3.0*v1-7.0/6.0*v2+11.0/6.0*v3  # r=2
    q1 = -1.0/6.0*v2+5.0/6.0*v3+1.0/3.0*v4  # r=1
    q2 = 1.0/3.0*v3+5.0/6.0*v4-1.0/6.0*v5  # r=0
    return w0 * q0 + w1 * q1 + w2 * q2

def _reconstruct_weno5_right(self, v1, v2, v3, v4, v5):
    eps = 1e-6
    beta0 = (13.0/12.0)*(v1 - 2*v2 + v3)**2 + (1.0/4.0)*(v1 - 4*v2 + 3*v3)**2
    beta1 = (13.0/12.0)*(v2 - 2*v3 + v4)**2 + (1.0/4.0)*(v2 - v4)**2
    beta2 = (13.0/12.0)*(v3 - 2*v4 + v5)**2 + (1.0/4.0)*(3*v3 - 4*v4 + v5)**2
    d0 = 3.0/10.0
    d1 = 3.0/5.0
    d2 = 1.0/10.0
    alpha0 = d0 / (eps + beta0)**2
    alpha1 = d1 / (eps + beta1)**2
    alpha2 = d2 / (eps + beta2)**2
    alpha = alpha0 + alpha1 + alpha2
    w0 = alpha0 / alpha
    w1 = alpha1 / alpha
    w2 = alpha2 / alpha
    q0 = -1.0/6.0*v1+5.0/6.0*v2+1.0/3.0*v3  # r=2
    q1 = 1.0/3.0*v2+5.0/6.0*v3-1.0/6.0*v4  # r=1
    q2 = 11.0/6.0*v3-7.0/6.0*v4+1.0/3.0*v5  # r=0
    return w0 * q0 + w1 * q1 + w2 * q2
