# config.py
class CfdConfig:
    def __init__(self):
        self.ic_type = "step"
        #self.ic_type = "sin"
        self.recon_scheme = "eno"  # 0=ENO, 1=WENO
        self.flux_type = "rusanov"  # "rusanov", "engquist-osher"
        #self.flux_type = "engquist-osher"  # "rusanov", "engquist-osher"
        self.rk_order = 1
        self.wave_speed = 1.0
        self.final_time = 0.625
        self.dt = 0.025
        
        self.boundary_type = "periodic"
        self.left_boundary_value = 1.0  # Dirichlet左边界值
        self.right_boundary_value = 2.0 # Dirichlet右边界值
        
        self.spatial_order = 2
        
    def with_reconstruction(self, scheme, order=None):
        """专用配置：重建方案（链式调用）"""
        self.recon_scheme = scheme.lower()  # 统一小写，避免大小写问题
        
        # 智能默认阶数
        if order is not None:
            self.spatial_order = order
        else:
            if self.recon_scheme.startswith("weno"):
                self.spatial_order = 5
            elif self.recon_scheme == "eno":
                self.spatial_order = 3  # ENO默认3阶
            else:
                raise ValueError(f"不支持的重建格式：{scheme}（仅支持 eno/weno）")
        
        return self
        
    def with_boundary(self, bc_type, left_value=None, right_value=None):
        self.boundary_type = bc_type
        if left_value is not None:
            self.left_boundary_value = left_value
        if right_value is not None:
            self.right_boundary_value = right_value
        return self   