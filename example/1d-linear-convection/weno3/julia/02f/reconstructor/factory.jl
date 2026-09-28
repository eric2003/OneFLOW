# julia/reconstructor/factory.jl

"""
ReconstructorFactory
对标 Python: ReconstructorFactory.create(config, domain)
"""
module ReconstructorFactory

# ✅ 不要用 using ..EnoReconstructor（它们不是模块）
# 直接通过 Main. 引用顶层定义的类型

function create(config::Any, domain::Any)
    scheme = lowercase(string(getfield_safe(config, :recon_scheme, "")))
    if scheme == ""
        error("必须先通过 with_reconstruction 设置重建格式！")
    end

    # 处理 WENO 默认命名（如 Python）
    if scheme == "weno"
        order = getfield_safe(config, :spatial_order, 5)
        scheme = "weno$(order)"
    end

    if scheme == "eno"
        order = getfield_safe(config, :spatial_order, 3)
        return Main.EnoReconstructor(order, domain.ntcells)
    elseif scheme == "weno3"
        return Main.Weno3Reconstructor()
    else
        error("不支持的重建格式: $scheme（仅支持 eno/weno3）")
    end
end

# 辅助函数（复制自 domain.jl）
function getfield_safe(obj, name::Symbol, default)
    return isdefined(obj, name) ? getfield(obj, name) : default
end

end  # module ReconstructorFactory

export ReconstructorFactory