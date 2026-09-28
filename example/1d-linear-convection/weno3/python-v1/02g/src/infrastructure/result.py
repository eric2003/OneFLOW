# src/infrastructure/result.py

class ResultAssembler:
    @staticmethod
    def assemble(cfd: 'Cfd') -> dict:
        u_numerical = cfd.solution.u[cfd.domain.ist:cfd.domain.ied].copy()
        
        # ✅ 从 problem 获取解析解（唯一正确路径）
        try:
            analytical = cfd.problem.exact_solution(cfd)
        except NotImplementedError:
            analytical = None  # 或 np.full_like(u_numerical, np.nan)

        return {
            "x": cfd.domain.mesh.xcc,
            "numerical": u_numerical,
            "analytical": analytical,
            "config": {
                "scheme": cfd.config.recon_scheme,
                "order": cfd.config.spatial_order,
                "rk_order": cfd.config.rk_order,
                "final_time": cfd.config.final_time,
                "problem": cfd.problem.name,
            }
        }