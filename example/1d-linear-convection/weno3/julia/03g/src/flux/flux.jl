# src/flux/flux.jl

# 加载组件（顺序很重要！）
include("base.jl")
include("rusanov.jl")
include("engquist_osher.jl")
include("factory.jl")

# 导出（如果你未来用模块）
# export InviscidFluxCalculator, RusanovFluxCalculator, EngquistOsherFluxCalculator, FluxCalculatorFactory