# src/reconstructor/factory.jl

module ReconstructorFactory

const _REGISTRY = Dict{String, Function}()

function register(name::String, ctor::Function)
    if haskey(_REGISTRY, name)
        @warn "重构器 '$name' 已注册，将被覆盖"
    end
    _REGISTRY[name] = ctor
end

function create(config::Any, domain::Any)
    scheme = lowercase(string(getfield_safe(config, :recon_scheme, "")))
    if scheme == ""
        error("必须先通过 with_reconstruction 设置重建格式！")
    end

    # 处理 WENO 默认命名
    if scheme == "weno"
        order = getfield_safe(config, :spatial_order, 5)
        scheme = "weno$(order)"
    end

    if !haskey(_REGISTRY, scheme)
        available = sort(collect(keys(_REGISTRY)))
        error("未注册的重构器: '$scheme'。可用选项: $available")
    end

    return _REGISTRY[scheme](config, domain)
end

# ✅ 手动注册（与你原有风格一致）
register("eno", (config, domain) -> Main.EnoReconstructor(config.spatial_order, domain.ntcells))
register("weno3", (config, domain) -> Main.Weno3Reconstructor())
register("weno5", (config, domain) -> Main.Weno5Reconstructor())  # ← 新增

# ---------------------- 辅助函数 ----------------------
function getfield_safe(obj, name::Symbol, default)
    return isdefined(obj, name) ? getfield(obj, name) : default
end

end  # module ReconstructorFactory

export ReconstructorFactory