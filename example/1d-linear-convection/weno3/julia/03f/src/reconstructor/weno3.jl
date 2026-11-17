# julia/reconstructor/weno3.jl
"""
WENO3 重构器（与 reconstructor/weno3.py 完全同构）
"""

mutable struct Weno3Reconstructor
    # 无字段，与 Python 一致
end

function reconstruct(rec::Weno3Reconstructor, q::Vector{Float64}, cfd::Any)
    domain = cfd.domain
    solution = cfd.solution
    _reconstruct_left_interfaces_weno3(domain, q, solution.q_face_left)
    _reconstruct_right_interfaces_weno3(domain, q, solution.q_face_right)
end

function _reconstruct_left_interfaces_weno3(domain, u, qL)
    """在每个 i+1/2 界面，计算左单元贡献的 qL (即 u_{i+1/2}^-)"""
    for i in (domain.ist - 1):(domain.ied - 1)
        j = i - (domain.ist - 1) + 1  # ← Julia 1-based: j = i - (ist-1) 对应 Python j = i - (ist-1)
        v1, v2, v3 = u[i-1], u[i], u[i+1]
        qL[j] = _reconstruct_from_left_biased_stencil(v1, v2, v3)
    end
end

function _reconstruct_right_interfaces_weno3(domain, u, qR)
    """在每个 i+1/2 界面，计算右单元贡献的 qR (即 u_{i+1/2}^+)"""
    for i in domain.ist:domain.ied
        j = i - domain.ist + 1  # ← Julia 1-based: j = i - ist 对应 Python j = i - ist
        v1, v2, v3 = u[i-1], u[i], u[i+1]
        qR[j] = _reconstruct_from_right_biased_stencil(v1, v2, v3)
    end
end

function _reconstruct_from_left_biased_stencil(v1, v2, v3)
    eps = 1e-6
    beta0 = (v2 - v1)^2
    beta1 = (v3 - v2)^2
    d0 = 1.0/3.0
    d1 = 2.0/3.0
    alpha0 = d0 / (eps + beta0)^2
    alpha1 = d1 / (eps + beta1)^2
    alpha = alpha0 + alpha1
    w0 = alpha0 / alpha
    w1 = alpha1 / alpha
    q0 = -0.5*v1 + 1.5*v2  # r=1 stencil
    q1 = 0.5*v2 + 0.5*v3   # r=0 stencil
    return w0 * q0 + w1 * q1
end

function _reconstruct_from_right_biased_stencil(v1, v2, v3)
    eps = 1e-6
    beta0 = (v2 - v1)^2
    beta1 = (v3 - v2)^2
    d0 = 2.0/3.0
    d1 = 1.0/3.0
    alpha0 = d0 / (eps + beta0)^2
    alpha1 = d1 / (eps + beta1)^2
    alpha = alpha0 + alpha1
    w0 = alpha0 / alpha
    w1 = alpha1 / alpha
    q0 = 0.5*v1 + 0.5*v2  # r=1 stencil
    q1 = 1.5*v2 - 0.5*v3  # r=0 stencil
    return w0 * q0 + w1 * q1
end