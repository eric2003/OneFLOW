# julia/solver.jl
"""
CFD 求解器主类（与 solver.py 完全同构）
"""

include("mesh.jl")
include("domain.jl")
include("solution.jl")
include("initial_condition.jl")
include("boundary.jl")
include("flux/flux.jl")
include("residual.jl")
include("time_integration.jl")
include("reconstructor/eno.jl")
include("reconstructor/weno3.jl")
include("reconstructor/factory.jl")

# 导入工厂模块（必须在顶层！）
using .FluxCalculatorFactory
using .TimeIntegratorFactory
using .BoundaryConditionFactory
using .ReconstructorFactory

# ---------------------- Cfd 求解器 ----------------------
mutable struct Cfd
    config::Any
    domain::Domain
    solution::Solution
    reconstructor::Any
    residual_calculator::ResidualCalculator
    integrator::Any
    boundary_condition::Any
    result::Dict{String, Any}
    
	function Cfd(config::Any, mesh::Mesh)
		domain = Domain(config, mesh)
		solution = Solution(config, domain)
		reconstructor = ReconstructorFactory.create(config, domain)
		
		# 1. 初始上下文（仅包含不依赖其他组件的字段）
		full_cfd = (
			config = config,
			domain = domain,
			solution = solution,
			reconstructor = reconstructor
			# 注意：不预先占位 nothing！
		)
		
		# 2. 创建 boundary_condition（只依赖 config + domain）
		boundary_condition = BoundaryConditionFactory.create(full_cfd)
		full_cfd = merge(full_cfd, (boundary_condition = boundary_condition,))
		
		# 3. 创建 residual_calculator（依赖上面所有字段）
		residual_calculator = ResidualCalculator(full_cfd)
		full_cfd = merge(full_cfd, (residual_calculator = residual_calculator,))
		
		# 4. 创建 integrator（依赖 residual_calculator 等）
		integrator = TimeIntegratorFactory.create(full_cfd)
		full_cfd = merge(full_cfd, (integrator = integrator,))
		
		# 5. 注入完整 self 到组件（确保它们能访问彼此）
		residual_calculator.cfd = full_cfd
		integrator.base.cfd = full_cfd
		
		result = Dict{String, Any}()
		new(config, domain, solution, reconstructor, residual_calculator, integrator, boundary_condition, result)
	end
end

"""
通用对流问题的解析解：u(x, T) = u0(x - c*T)，周期边界
"""
function exact_solution(cfd::Cfd)
    x = cfd.domain.mesh.xcc
    T = cfd.config.final_time
    c = cfd.config.wave_speed
    L = cfd.domain.mesh.L
    
    # 周期平移：确保在 [x0, x0 + L) 内
    x_shifted = @. (x - c * T + L) % L
    
    # ✅ 使用工厂创建初始条件（对标 Python）
    ic = InitialConditionFactory.create(cfd.config.ic_type, cfd.config)
    
    return evaluate_at(ic, x_shifted)
end

"""
主求解循环
"""
function run!(cfd::Cfd)
    # 应用初始边界条件并同步 old field
    apply!(cfd.boundary_condition, cfd.solution.u)
    update_old_field(cfd.solution)
    
    t = 0.0
    dt_old = cfd.config.dt
    dt = dt_old
    
    while t < cfd.config.final_time
        if t + dt > cfd.config.final_time
            dt = cfd.config.final_time - t
        end
		#@show t, dt, maximum(cfd.solution.u), minimum(cfd.solution.u)
        # 执行时间步
        step(cfd.integrator, dt)
        t += dt
    end
    
    # 恢复 dt
    cfd.config.dt = dt_old
    
    # 整理结果
    u_numerical = cfd.solution.u[cfd.domain.ist:cfd.domain.ied-1]  # Python: [ist:ied]
    analytical = exact_solution(cfd)
    
    cfd.result = Dict(
        "x" => cfd.domain.mesh.xcc,
        "numerical" => u_numerical,
        "analytical" => analytical,
        "config" => Dict(
            "scheme" => cfd.config.recon_scheme,
            "order" => cfd.config.spatial_order,
            "rk_order" => cfd.config.rk_order,
            "final_time" => cfd.config.final_time
        )
    )
    
    return u_numerical
end