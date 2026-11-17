from abc import ABC, abstractmethod

# ---------------------- 1. 抽象通量计算基类（统一接口） ----------------------
class InviscidFluxCalculator(ABC):
    """无粘通量计算抽象基类：定义一维CFD通量计算接口"""
    def __init__(self, cfd):
        self.cfd = cfd
        self.config = cfd.config
        self.mesh = cfd.domain.mesh
        self.wave_speed = self.config.wave_speed

    @abstractmethod
    def compute(self, q_face_left, q_face_right, flux):
        """
        计算无粘通量（核心接口）
        :param q_face_left: 左界面值数组
        :param q_face_right: 右界面值数组
        :param flux: 输出通量数组
        :return: None
        """
        pass

# ---------------------- 2. 具体通量计算子类（隔离不同格式） ----------------------
class RusanovFluxCalculator(InviscidFluxCalculator):
    """Rusanov（Lax-Friedrichs）通量"""
    def compute(self, q_face_left, q_face_right, flux):
        for i in range(self.mesh.nnodes):
            u_L = q_face_left[i]
            u_R = q_face_right[i]
            c_L = self.wave_speed
            c_R = self.wave_speed
            F_L = c_L * u_L  # Flux from left state
            F_R = c_R * u_R  # Flux from right state
            Smax = max(abs(c_L), abs(c_R))  # Maximum wave speed
            flux[i] = 0.5 * (F_L + F_R) - 0.5 * Smax * (u_R - u_L)

class EngquistOsherFluxCalculator(InviscidFluxCalculator):
    """Engquist-Osher通量（线性对流专用）"""
    def compute(self, q_face_left, q_face_right, flux):
        for i in range(self.mesh.nnodes):
            c = self.wave_speed
            cp = 0.5 * (c + abs(c))
            cm = 0.5 * (c - abs(c))
            u_L = q_face_left[i]
            u_R = q_face_right[i]
            flux[i] = cp * u_L + cm * u_R

# ---------------------- 3. 通量计算器工厂（统一创建逻辑） ----------------------
class FluxCalculatorFactory:
    @staticmethod
    def create(cfd):
        """根据配置创建通量计算器实例"""
        flux_type = cfd.config.flux_type
        flux_mapping = {
            0: RusanovFluxCalculator,
            1: EngquistOsherFluxCalculator,
            # 新增通量格式只需加键值对：2: LaxWendroffFluxCalculator
        }
        if flux_type not in flux_mapping:
            raise ValueError(f"不支持的通量类型：{flux_type}（可选：{list(flux_mapping.keys())}）")
        return flux_mapping[flux_type](cfd)