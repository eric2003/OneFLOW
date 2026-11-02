# src/flux/factory.jl

module FluxCalculatorFactory

# 全局注册表
const _REGISTRY = Dict{String, Function}()

function register(name::String, ctor::Function)
    if haskey(_REGISTRY, name)
        @warn "通量计算器 '$name' 已注册，将被覆盖"
    end
    _REGISTRY[name] = ctor
end

function create(cfd::Any)
    flux_type = cfd.config.flux_type
    if !haskey(_REGISTRY, flux_type)
        error("未注册的通量计算器: '$flux_type'。可用: $(keys(_REGISTRY))")
    end
    return _REGISTRY[flux_type](cfd)
end

register("rusanov", cfd -> Main.RusanovFluxCalculator(cfd))
register("engquist-osher", cfd -> Main.EngquistOsherFluxCalculator(cfd))

end # module FluxCalculatorFactory