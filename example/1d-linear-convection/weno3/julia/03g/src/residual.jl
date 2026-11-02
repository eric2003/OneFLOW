# julia/residual.jl
"""
残差计算器（与 residual.py 完全同构）
- 封装重建→通量→散度完整流程
- 依赖 cfd 的多个字段
"""

include("mesh.jl")

mutable struct ResidualCalculator
    cfd::Any
    config::Any
    domain::Any
    solution::Any
    mesh::Mesh
    reconstructor::Any
    flux_calculator::Any
    
	function ResidualCalculator(cfd::Any)
		config = cfd.config
		domain = cfd.domain
		solution = cfd.solution
		mesh = domain.mesh
		reconstructor = cfd.reconstructor
		
		# ✅ 内部创建 flux_calculator（对标 Python）
		flux_calculator = FluxCalculatorFactory.create(cfd)
		
		new(cfd, config, domain, solution, mesh, reconstructor, flux_calculator)
	end
end


"""
计算完整残差（对外唯一接口）
"""
function compute!(calc::ResidualCalculator)
    _reconstruct(calc)
    _compute_inviscid_flux(calc)
    _compute_flux_divergence(calc)
end

"""
私有方法：界面值重建
"""
function _reconstruct(calc::ResidualCalculator)
    reconstruct(calc.reconstructor, calc.solution.u, calc.cfd)
end

"""
私有方法：计算无粘通量
"""
function _compute_inviscid_flux(calc::ResidualCalculator)
    compute!(calc.flux_calculator,
             calc.solution.q_face_left,
             calc.solution.q_face_right,
             calc.solution.flux)
end

"""
私有方法：计算通量散度（残差 = -dF/dx）
"""
function _compute_flux_divergence(calc::ResidualCalculator)
    solution = calc.solution
    mesh = calc.mesh
    for i in 1:mesh.ncells
        solution.res[i] = -(solution.flux[i+1] - solution.flux[i]) / mesh.dx
    end
end