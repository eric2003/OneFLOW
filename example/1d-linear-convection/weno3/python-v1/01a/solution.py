# solution.py
import numpy as np
from initial_condition import InitialConditionFactory

class Solution:
    def __init__(self, config, domain):
        """
        初始化求解过程中的动态数据
        :param domain: Domain实例（用于确定数组尺寸）
        """
        self.domain = domain
        mesh = domain.mesh

        # 界面值和通量（维度依赖mesh.nnodes）
        self.q_face_left = np.zeros(mesh.nnodes)   # 左界面值
        self.q_face_right = np.zeros(mesh.nnodes)  # 右界面值
        self.flux = np.zeros(mesh.nnodes)          # 通量

        # 残差（维度依赖mesh.ncells）
        self.res = np.zeros(mesh.ncells)           # 残差

        # 解数组（维度依赖ntcells，含ghost层）
        self.u = np.zeros(domain.ntcells)          # 当前解
        self.un = np.zeros(domain.ntcells)         # 上一时间步解
        
    def reset_solution(self):
        """重置解数组为初始状态"""
        self.u.fill(0.0)
        self.un.fill(0.0)
        
    def initialize_from_config(self, config):
        """根据配置初始化场"""
        ic = InitialConditionFactory.create(config)
        ic.apply(self)
        
    def update_old_field(self):
        """更新旧场"""
        self.un[:] = self.u[:]