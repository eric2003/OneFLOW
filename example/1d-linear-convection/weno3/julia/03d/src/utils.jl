# src/utils.jl

# ---------------------- 辅助函数 ----------------------
"""
模拟 Python 的 getattr(obj, name, default)
"""
function getfield_safe(obj, name::Symbol, default)
    return isdefined(obj, name) ? getfield(obj, name) : default
end
