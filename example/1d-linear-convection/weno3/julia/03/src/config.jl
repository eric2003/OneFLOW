# julia/config.jl
"""
CfdConfig：与 Python config.py 完全同构
"""
mutable struct CfdConfig
    ic_type::String
    recon_scheme::String
    flux_type::String
    rk_order::Int
    wave_speed::Float64
    final_time::Float64
    dt::Float64
    boundary_type::String
    left_boundary_value::Float64
    right_boundary_value::Float64
    spatial_order::Int
    
    function CfdConfig()
        new(
            "step",
            "eno",
            "rusanov",
            1,
            1.0,
            0.625,
            0.025,
            "periodic",
            1.0,
            2.0,
            2
        )
    end
end

"""
专用配置：重建方案（链式调用）
"""
function with_reconstruction(cfg::CfdConfig, scheme::String, order::Union{Int, Nothing}=nothing)
    cfg.recon_scheme = lowercase(scheme)
    
    if order !== nothing
        cfg.spatial_order = order
    else
        if startswith(cfg.recon_scheme, "weno")
            cfg.spatial_order = 5
        elseif cfg.recon_scheme == "eno"
            cfg.spatial_order = 3
        else
            error("不支持的重建格式：$scheme（仅支持 eno/weno）")
        end
    end
    
    return cfg  # 支持链式调用
end

"""
专用配置：边界条件（链式调用）
"""
function with_boundary(cfg::CfdConfig, bc_type::String; left_value=nothing, right_value=nothing)
    cfg.boundary_type = bc_type
    if left_value !== nothing
        cfg.left_boundary_value = left_value
    end
    if right_value !== nothing
        cfg.right_boundary_value = right_value
    end
    return cfg
end

# ---------------------- 辅助函数 ----------------------
"""
模拟 Python 的 getattr(obj, name, default)
"""
function getfield_safe(obj, name::Symbol, default)
    return isdefined(obj, name) ? getfield(obj, name) : default
end