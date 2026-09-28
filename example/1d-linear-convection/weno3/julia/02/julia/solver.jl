# julia/solver.jl
"""
CFD 求解器主类（与 solver.py 完全同构）
"""

include("mesh.jl")
include("domain.jl")
include("solution.jl")
include("initial_condition.jl")
include("boundary.jl")
include("flux.jl")
include("residual.jl")
include("time_integration.jl")
include("reconstructor/eno.jl")
include("reconstructor/weno3.jl")

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
		
		# 在 Cfd 构造函数中
		# =============== 重建器创建 ===============
		recon_scheme = config.recon_scheme
		spatial_order = config.spatial_order

		# ✅ 模拟 Python ReconstructorFactory 的逻辑
		if recon_scheme == "weno"
			recon_scheme = "weno$(spatial_order)"
		end

		# 现在根据 recon_scheme 创建
		reconstructor = if recon_scheme == "eno"
			EnoReconstructor(spatial_order, domain.ntcells)
		elseif recon_scheme == "weno3"
			Weno3Reconstructor()
		else
			error("不支持的重建格式: $recon_scheme")
		end
        
        # 通量计算器
        flux_calculator = if config.flux_type == "rusanov"
            RusanovFluxCalculator((config=config, domain=domain))
        elseif config.flux_type == "engquist-osher"
            EngquistOsherFluxCalculator((config=config, domain=domain))
        else
            error("不支持的通量类型: $(config.flux_type)")
        end
        
        # 残差计算器
        residual_calculator = ResidualCalculator(
            (config=config, domain=domain, solution=solution, reconstructor=reconstructor),
            flux_calculator
        )
        
        # 边界条件
        boundary_condition = if config.boundary_type == "periodic"
            PeriodicBoundary((config=config, domain=domain))
        elseif config.boundary_type == "dirichlet"
            DirichletBoundary((config=config, domain=domain))
        elseif config.boundary_type == "neumann"
            NeumannBoundary((config=config, domain=domain))
        else
            error("不支持的边界类型: $(config.boundary_type)")
        end
        
		# 构造用于积分器初始化的上下文对象（NamedTuple，模拟 cfd 接口）
		integrator_context = (
			config = config,
			domain = domain,
			solution = solution,
			residual_calculator = residual_calculator,
			boundary_condition = boundary_condition
		)

		# 使用注册表工厂创建时间推进器
		integrator = create_integrator(integrator_context)		
		
		#@show typeof(integrator)
        
        # 注入 cfd 到 residual_calculator 和 integrator
        residual_calculator.cfd = (config=config, domain=domain, solution=solution, reconstructor=reconstructor, residual_calculator=residual_calculator, integrator=integrator, boundary_condition=boundary_condition)
        integrator.base.cfd = residual_calculator.cfd
        
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
    
    # 获取 IC 实例并评估
    ic = if cfd.config.ic_type == "step"
        StepFunctionIC(cfd.config)
    elseif cfd.config.ic_type == "sin"
        SineWaveIC(cfd.config)
    elseif cfd.config.ic_type == "gaussian"
        GaussianPulseIC(cfd.config)
    else
        error("未知初始条件: $(cfd.config.ic_type)")
    end
    
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