# julia/reconstructor/eno.jl
"""
ENO 重构器（与 reconstructor/eno.py 完全同构）
"""

# ---------------------- ENO 系数初始化 ----------------------
function _init_eno_coef!(spatial_order::Int, coef::Matrix{Float64})
    if spatial_order == 1
        coef[1, 1] = 1.0
        coef[2, 1] = 1.0
    elseif spatial_order == 2
        coef[1, 1:2] = [3.0/2.0, -1.0/2.0]
        coef[2, 1:2] = [1.0/2.0,  1.0/2.0]
        coef[3, 1:2] = [-1.0/2.0, 3.0/2.0]
    elseif spatial_order == 3
        coef[1, 1:3] = [11.0/6.0, -7.0/6.0, 1.0/3.0]
        coef[2, 1:3] = [1.0/3.0,  5.0/6.0, -1.0/6.0]
        coef[3, 1:3] = [-1.0/6.0, 5.0/6.0, 1.0/3.0]
        coef[4, 1:3] = [1.0/3.0, -7.0/6.0, 11.0/6.0]
    elseif spatial_order == 4
        coef[1, 1:4] = [25.0/12.0, -23.0/12.0, 13.0/12.0, -1.0/4.0]
        coef[2, 1:4] = [1.0/4.0, 13.0/12.0, -5.0/12.0, 1.0/12.0]
        coef[3, 1:4] = [-1.0/12.0, 7.0/12.0, 7.0/12.0, -1.0/12.0]
        coef[4, 1:4] = [1.0/12.0, -5.0/12.0, 13.0/12.0, 1.0/4.0]
        coef[5, 1:4] = [-1.0/4.0, 13.0/12.0, -23.0/12.0, 25.0/12.0]
    elseif spatial_order == 5
        coef[1, 1:5] = [137.0/60.0, -163.0/60.0, 137.0/60.0, -21.0/20.0, 1.0/5.0]
        coef[2, 1:5] = [1.0/5.0, 77.0/60.0, -43.0/60.0, 17.0/60.0, -1.0/20.0]
        coef[3, 1:5] = [-1.0/20.0, 9.0/20.0, 47.0/60.0, -13.0/60.0, 1.0/30.0]
        coef[4, 1:5] = [1.0/30.0, -13.0/60.0, 47.0/60.0, 9.0/20.0, -1.0/20.0]
        coef[5, 1:5] = [-1.0/20.0, 17.0/60.0, -43.0/60.0, 77.0/60.0, 1.0/5.0]
        coef[6, 1:5] = [1.0/5.0, -21.0/20.0, 137.0/60.0, -163.0/60.0, 137.0/60.0]
    elseif spatial_order == 6
        coef[1, 1:6] = [49.0/20.0, -71.0/20.0, 79.0/20.0, -163.0/60.0, 31.0/30.0, -1.0/6.0]
        coef[2, 1:6] = [1.0/6.0, 29.0/20.0, -21.0/20.0, 37.0/60.0, -13.0/60.0, 1.0/30.0]
        coef[3, 1:6] = [-1.0/30.0, 11.0/30.0, 19.0/20.0, -23.0/60.0, 7.0/60.0, -1.0/60.0]
        coef[4, 1:6] = [1.0/60.0, -2.0/15.0, 37.0/60.0, 37.0/60.0, -2.0/15.0, 1.0/60.0]
        coef[5, 1:6] = [-1.0/60.0, 7.0/60.0, -23.0/60.0, 19.0/20.0, 11.0/30.0, -1.0/30.0]
        coef[6, 1:6] = [1.0/30.0, -13.0/60.0, 37.0/60.0, -21.0/20.0, 29.0/20.0, 1.0/6.0]
        coef[7, 1:6] = [-1.0/6.0, 31.0/30.0, -163.0/60.0, 79.0/20.0, -71.0/20.0, 49.0/20.0]
    elseif spatial_order == 7
        coef[1, 1:7] = [363.0/140.0, -617.0/140.0, 853.0/140.0, -2341.0/420.0, 667.0/210.0, -43.0/42.0, 1.0/7.0]
        coef[2, 1:7] = [1.0/7.0, 223.0/140.0, -197.0/140.0, 153.0/140.0, -241.0/420.0, 37.0/210.0, -1.0/42.0]
        coef[3, 1:7] = [-1.0/42.0, 13.0/42.0, 153.0/140.0, -241.0/420.0, 109.0/420.0, -31.0/420.0, 1.0/105.0]
        coef[4, 1:7] = [1.0/105.0, -19.0/210.0, 107.0/210.0, 319.0/420.0, -101.0/420.0, 5.0/84.0, -1.0/140.0]
        coef[5, 1:7] = [-1.0/140.0, 5.0/84.0, -101.0/420.0, 319.0/420.0, 107.0/210.0, -19.0/210.0, 1.0/105.0]
        coef[6, 1:7] = [1.0/105.0, -31.0/420.0, 109.0/420.0, -241.0/420.0, 153.0/140.0, 13.0/42.0, -1.0/42.0]
        coef[7, 1:7] = [-1.0/42.0, 37.0/210.0, -241.0/420.0, 153.0/140.0, -197.0/140.0, 223.0/140.0, 1.0/7.0]
        coef[8, 1:7] = [1.0/7.0, -43.0/42.0, 667.0/210.0, -2341.0/420.0, 853.0/140.0, -617.0/140.0, 363.0/140.0]
    else
        error("ENO 系数未实现 order=$spatial_order")
    end
end

# ---------------------- ENO 重构器 ----------------------
mutable struct EnoReconstructor
    cfd::Any
    config::Any
    domain::Any
    spatial_order::Int
    ntcells::Int
    lmc::Vector{Int}
    coef::Matrix{Float64}
    dd::Matrix{Float64}
    
    function EnoReconstructor(cfd::Any)
        config = cfd.config
        domain = cfd.domain
        spatial_order = config.spatial_order
        ntcells = domain.ntcells
        lmc = zeros(Int, ntcells)
        coef = zeros(Float64, spatial_order + 1, spatial_order)
        dd = zeros(Float64, spatial_order, ntcells)
        _init_eno_coef!(spatial_order, coef)
        new(cfd, config, domain, spatial_order, ntcells, lmc, coef, dd)
    end
end

function reconstruct(rec::EnoReconstructor, q::Vector{Float64}, cfd::Any)
    # 1. 差商计算 (dd[1,:] = q)
    @views rec.dd[1, :] .= q
    for m in 2:rec.spatial_order
        for j in 1:(rec.ntcells - m + 1)
            rec.dd[m, j] = rec.dd[m-1, j+1] - rec.dd[m-1, j]
        end
    end

    # 2. 选择 smoothest stencil
    domain = cfd.domain
    for i in (domain.ist - 1):(domain.ied)  # Python: range(ist-1, ied+1) → ied+1-1 = ied
        rec.lmc[i] = i
        for m in 2:rec.spatial_order
            if abs(rec.dd[m, rec.lmc[i] - 1]) < abs(rec.dd[m, rec.lmc[i]])
                rec.lmc[i] -= 1
            end
        end
    end

    # 3. 重构界面值
    solution = cfd.solution
    for i in domain.ist:(domain.ied)  # Python: range(ist, ied+1) → ied+1-1 = ied
        j = i - domain.ist + 1  # Julia 1-based
        k1 = rec.lmc[i - 1]
        k2 = rec.lmc[i]
        r1 = (i - 1) - k1 + 1
        r2 = i - k2 + 1
        solution.q_face_left[j] = 0.0
        solution.q_face_right[j] = 0.0
        for m in 1:rec.spatial_order
            solution.q_face_left[j] += q[k1 + m - 1] * rec.coef[r1 + 1, m]
            solution.q_face_right[j] += q[k2 + m - 1] * rec.coef[r2, m]
        end
    end
end