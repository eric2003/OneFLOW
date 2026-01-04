# flux/base.py
"""
抽象通量计算基类
"""

from abc import ABC, abstractmethod

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