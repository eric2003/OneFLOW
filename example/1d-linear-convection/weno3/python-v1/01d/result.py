# result.py

from initial_condition import InitialConditionFactory

class ResultAssembler:
    @staticmethod
    def assemble(cfd: 'Cfd') -> dict:
        u_numerical = cfd.solution.u[cfd.domain.ist:cfd.domain.ied].copy()
        analytical = InitialConditionFactory.get_exact_solution(cfd)
        
        return {
            "x": cfd.domain.mesh.xcc,
            "numerical": u_numerical,
            "analytical": analytical,
            "config": {
                "scheme": cfd.config.recon_scheme,
                "order": cfd.config.spatial_order,
                "rk_order": cfd.config.rk_order,
                "final_time": cfd.config.final_time
            }
        }