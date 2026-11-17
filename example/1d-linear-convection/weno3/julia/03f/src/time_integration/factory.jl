# src/time_integration/factory.jl

module TimeIntegratorFactory

const _REGISTRY = Dict{String, Function}()

function register(name::String, ctor::Function)
    if haskey(_REGISTRY, name)
        @warn "时间积分器 '$name' 已注册，将被覆盖"
    end
    _REGISTRY[name] = ctor
end

function create(cfd::Any)
    rk_order = cfd.config.rk_order
    name = "rk$rk_order"
    if !haskey(_REGISTRY, name)
        available = sort(collect(keys(_REGISTRY)))
        error("未注册的时间积分器: '$name'。可用选项: $available")
    end
    return _REGISTRY[name](cfd)
end

# ✅ 关键：用 Main. 引用在 Main 模块中定义的类型
register("rk1", cfd -> Main.RK1Integrator(cfd))
register("rk2", cfd -> Main.RK2Integrator(cfd))
register("rk3", cfd -> Main.RK3Integrator(cfd))

end  # module TimeIntegratorFactory