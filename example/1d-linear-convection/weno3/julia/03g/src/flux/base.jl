# src/flux/base.jl
"""
抽象通量计算器接口
Julia 无 ABC，用文档约定：
- 所有子类型必须实现 `compute!(calc, qL, qR, flux)`
"""
abstract type InviscidFluxCalculator end