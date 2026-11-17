# julia/domain.jl
"""
Domain 结构体（与 domain.py 完全同构）
- 保存 config
- nghosts 计算逻辑 1:1
- ist/ied 为 0-based（与 Python 一致）
"""

include("mesh.jl")

mutable struct Domain
    config::Any      # 保存 config（与 Python 一致）
    mesh::Mesh
    nghosts::Int
    ist::Int         # 0-based
    ied::Int         # 0-based
    ntcells::Int
end

function Domain(config::Any, mesh::Mesh)
    nghosts = _calc_nghosts(config)
    ist = nghosts + 1          # ← 1-based 起始索引
    ied = ist + mesh.ncells    # ← 1-based 结束索引（不包含）
    ntcells = mesh.ncells + 2 * nghosts
    Domain(config, mesh, nghosts, ist, ied, ntcells)
end

function _calc_nghosts(config::Any)
    scheme = lowercase(string(getfield_safe(config, :recon_scheme, "")))
    order = getfield_safe(config, :spatial_order, 2)

    if scheme == ""
        error("必须先通过 with_reconstruction 设置重建格式！")
    end

    if scheme == "eno"
        nghosts = order
    elseif startswith(scheme, "weno")
        nghosts = order ÷ 2 + 1
    else
        error("未知重建格式 $(scheme)，无法计算ghost层！")
    end

    if nghosts <= 0
        error("计算得到的ghost层数量无效：$(nghosts)（阶数$(order)，格式$(scheme)）")
    end

    return nghosts
end

function is_physical_cell(domain::Domain, idx::Int)
    return domain.ist <= idx < domain.ied
end

function get_physical_indices(domain::Domain)
    return domain.ist:(domain.ied - 1)
end
