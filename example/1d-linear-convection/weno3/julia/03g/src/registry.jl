"""
统一注册系统 + 通用工厂
- ComponentRegistry: 组件注册表
- @register_component: 装饰器宏
- BaseFactory: 通用工厂接口
"""

module ComponentRegistry

# 注册表：Dict{category => Dict{name => constructor}}
const _REGISTRIES = Dict{String, Dict{String, Function}}()
const _VERBOSE = Ref(true)

# ---------------------- 注册核心逻辑 ----------------------
"""
    register(category::String, name::String, ctor::Function)

注册一个组件构造函数到指定类别。
如果已存在同名组件且不同，则发出警告。
"""
function register(category::String, name::String, ctor::Function)
    if !haskey(_REGISTRIES, category)
        _REGISTRIES[category] = Dict{String, Function}()
    end

    registry = _REGISTRIES[category]
    if haskey(registry, name)
        if registry[name] !== ctor && _VERBOSE[]
            @warn "覆盖注册: $category.$name"
        end
    end

    registry[name] = ctor
    if _VERBOSE[]
        println("✅ 已注册: $category.$name -> $(nameof(ctor))")
    end
end

"""
    create(category::String, name::String, args...; kwargs...)

从注册表创建组件实例。
"""
function create(category::String, name::String, args...; kwargs...)
    if !haskey(_REGISTRIES, category)
        error("❌ 未知类别: $category (可用: $(collect(keys(_REGISTRIES))))")
    end

    registry = _REGISTRIES[category]
    lname = lowercase(name)
    if !haskey(registry, lname)
        available = sort(collect(keys(registry)))
        error("❌ 未找到: $category.$name (可用: $available)")
    end

    return registry[lname](args...; kwargs...)
end

"""
    list_all()

返回所有已注册组件（按类别）。
"""
function list_all()
    return Dict(cat => sort(collect(keys(reg))) for (cat, reg) in _REGISTRIES)
end

"""
    set_verbose(flag::Bool)

开启/关闭注册提示。
"""
function set_verbose(flag::Bool)
    _VERBOSE[] = flag
end

end  # module ComponentRegistry


# ---------------------- 装饰器宏：@register_component ----------------------
"""
@register_component(category, [name])

用法：
  @register_component("boundary", "periodic")
  struct PeriodicBoundary ...

若省略 name，则使用类型名的小写形式。
"""
macro register_component(category::String, name_expr)
    error("@register_component 需在类型定义前使用，且必须在模块顶层")
end

macro register_component(category::String)
    error("@register_component(category, name) 需指定 name 或在类型后使用")
end

# 重载：@register_component("category", "name") struct X ... end
macro register_component(category::String, name::String, ex)
    if !Meta.isexpr(ex, :struct)
        error("@register_component 必须作用于 struct 定义")
    end

    struct_name = ex.args[2]
    if Meta.isexpr(struct_name, :curly)
        struct_name = struct_name.args[1]
    end

    # 插入注册调用（在模块顶层）
    quote
        $(esc(ex))
        $(ComponentRegistry).register($(category), $(name), $(esc(struct_name)))
    end
end

# 重载：@register_component("category") struct X ... end → name = lowercase(nameof(X))
macro register_component(category::String, ex)
    if !Meta.isexpr(ex, :struct)
        error("@register_component 必须作用于 struct 定义")
    end

    struct_name = ex.args[2]
    if Meta.isexpr(struct_name, :curly)
        struct_name = struct_name.args[1]
    end

    name_str = string(struct_name) |> lowercase

    quote
        $(esc(ex))
        $(ComponentRegistry).register($(category), $(name_str), $(esc(struct_name)))
    end
end


# ---------------------- 通用工厂 ----------------------
module BaseFactory

using ..ComponentRegistry

"""
    create_component(category::String, name::String, args...; kwargs...)

通用工厂接口，与 Python BaseFactory.create_component 行为一致。
"""
function create_component(category::String, name::String, args...; kwargs...)
    return ComponentRegistry.create(category, name, args...; kwargs...)
end

"""
    get_available_components(category::String)

列出某类别下所有可用组件。
"""
function get_available_components(category::String)
    all = ComponentRegistry.list_all()
    return get(all, category, String[])
end

end  # module BaseFactory


# 导出接口
export ComponentRegistry, BaseFactory, @register_component