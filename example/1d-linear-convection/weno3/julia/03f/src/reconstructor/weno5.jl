# src/reconstructor/weno5.jl
"""
WENO5 重构器（与提供的 Python weno5.py 完全同构）
- 无字段（与 Weno3Reconstructor 一致）
- 所有逻辑 1:1 对齐 Python
"""

mutable struct Weno5Reconstructor
    # 无字段，与 Python class Weno5Reconstructor 一致
end

function reconstruct(rec::Weno5Reconstructor, q::Vector{Float64}, cfd::Any)
    domain = cfd.domain
    solution = cfd.solution
    _reconstruct_left_interfaces_weno5(domain, q, solution.q_face_left)
    _reconstruct_right_interfaces_weno5(domain, q, solution.q_face_right)
end

function _reconstruct_left_interfaces_weno5(domain, u, qL)
    """在每个 i+1/2 界面，计算左单元贡献的 qL (即 u_{i+1/2}^-)"""
    for i in (domain.ist - 1):(domain.ied - 1)
        j = i - (domain.ist - 1) + 1  # Julia 1-based
        v1, v2, v3, v4, v5 = u[i-2], u[i-1], u[i], u[i+1], u[i+2]
        qL[j] = _reconstruct_from_left_biased_stencil(v1, v2, v3, v4, v5)
    end
end

function _reconstruct_right_interfaces_weno5(domain, u, qR)
    """在每个 i+1/2 界面，计算右单元贡献的 qR (即 u_{i+1/2}^+)"""
    for i in domain.ist:domain.ied
        j = i - domain.ist + 1  # Julia 1-based
        v1, v2, v3, v4, v5 = u[i-2], u[i-1], u[i], u[i+1], u[i+2]
        qR[j] = _reconstruct_from_right_biased_stencil(v1, v2, v3, v4, v5)
    end
end

function _reconstruct_from_left_biased_stencil(v1, v2, v3, v4, v5)
    eps = 1e-6
    beta0 = (13.0/12.0)*(v1 - 2*v2 + v3)^2 + (1.0/4.0)*(v1 - 4*v2 + 3*v3)^2
    beta1 = (13.0/12.0)*(v2 - 2*v3 + v4)^2 + (1.0/4.0)*(v2 - v4)^2
    beta2 = (13.0/12.0)*(v3 - 2*v4 + v5)^2 + (1.0/4.0)*(3*v3 - 4*v4 + v5)^2
    d0 = 1.0/10.0
    d1 = 3.0/5.0
    d2 = 3.0/10.0
    alpha0 = d0 / (eps + beta0)^2
    alpha1 = d1 / (eps + beta1)^2
    alpha2 = d2 / (eps + beta2)^2
    alpha = alpha0 + alpha1 + alpha2
    w0 = alpha0 / alpha
    w1 = alpha1 / alpha
    w2 = alpha2 / alpha
    q0 = 1.0/3.0*v1 - 7.0/6.0*v2 + 11.0/6.0*v3  # r=2
    q1 = -1.0/6.0*v2 + 5.0/6.0*v3 + 1.0/3.0*v4  # r=1
    q2 = 1.0/3.0*v3 + 5.0/6.0*v4 - 1.0/6.0*v5  # r=0
    return w0 * q0 + w1 * q1 + w2 * q2
end

function _reconstruct_from_right_biased_stencil(v1, v2, v3, v4, v5)
    eps = 1e-6
    beta0 = (13.0/12.0)*(v1 - 2*v2 + v3)^2 + (1.0/4.0)*(v1 - 4*v2 + 3*v3)^2
    beta1 = (13.0/12.0)*(v2 - 2*v3 + v4)^2 + (1.0/4.0)*(v2 - v4)^2
    beta2 = (13.0/12.0)*(v3 - 2*v4 + v5)^2 + (1.0/4.0)*(3*v3 - 4*v4 + v5)^2
    d0 = 3.0/10.0
    d1 = 3.0/5.0
    d2 = 1.0/10.0
    alpha0 = d0 / (eps + beta0)^2
    alpha1 = d1 / (eps + beta1)^2
    alpha2 = d2 / (eps + beta2)^2
    alpha = alpha0 + alpha1 + alpha2
    w0 = alpha0 / alpha
    w1 = alpha1 / alpha
    w2 = alpha2 / alpha
    q0 = -1.0/6.0*v1 + 5.0/6.0*v2 + 1.0/3.0*v3  # r=2
    q1 = 1.0/3.0*v2 + 5.0/6.0*v3 - 1.0/6.0*v4  # r=1
    q2 = 11.0/6.0*v3 - 7.0/6.0*v4 + 1.0/3.0*v5  # r=0
    return w0 * q0 + w1 * q1 + w2 * q2
end