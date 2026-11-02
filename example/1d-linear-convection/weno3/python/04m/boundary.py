from abc import ABC, abstractmethod

# ---------------------- 边界条件抽象基类（统一接口） ----------------------
class BoundaryCondition(ABC):
    """边界条件抽象基类：定义所有边界条件必须实现的接口"""
    def __init__(self, cfd):
        self.cfd = cfd
        self.domain = cfd.domain
        self.config = cfd.config  # 可从配置读取边界参数（如进口速度、固壁温度等）
    
    @abstractmethod
    def apply(self, u):
        """
        应用边界条件到解数组
        :param u: 包含ghost层的解数组（会直接修改该数组）
        :return: None
        """
        pass

# ---------------------- 具体边界条件实现（可无限扩展） ----------------------
class PeriodicBoundary(BoundaryCondition):
    """周期边界条件（1D专用）"""
    def apply(self, u):
        nghosts = self.domain.nghosts
        ist = self.domain.ist
        ied = self.domain.ied
        
        # 左ghost层 = 右物理层
        for ig in range(nghosts):
            u[ist - 1 - ig] = u[ied - 1 - ig]
        
        # 右ghost层 = 左物理层
        for ig in range(nghosts):
            u[ied + ig] = u[ist + ig]

class DirichletBoundary(BoundaryCondition):
    """Dirichlet（固定值）边界条件（如进口固定速度、固壁零速度）"""
    def apply(self, u):
        nghosts = self.domain.nghosts
        ist = self.domain.ist
        ied = self.domain.ied
        
        # 左边界（进口）固定值（从配置读取）
        left_value = self.config.get("left_boundary_value", 1.0)
        for ig in range(nghosts):
            u[ist - 1 - ig] = left_value
        
        # 右边界（出口）固定值（从配置读取）
        right_value = self.config.get("right_boundary_value", 2.0)
        for ig in range(nghosts):
            u[ied + ig] = right_value

class NeumannBoundary(BoundaryCondition):
    """Neumann（零梯度）边界条件（如出口无梯度）"""
    def apply(self, u):
        nghosts = self.domain.nghosts
        ist = self.domain.ist
        ied = self.domain.ied
        
        # 左边界零梯度
        for ig in range(nghosts):
            u[ist - 1 - ig] = u[ist + ig]
        
        # 右边界零梯度
        for ig in range(nghosts):
            u[ied + ig] = u[ied - 1 - ig]

# ---------------------- 边界条件工厂（动态创建实例） ----------------------
class BoundaryConditionFactory:
    """边界条件工厂：根据配置创建对应边界条件实例"""
    @staticmethod
    def create(cfd):
        # 从配置读取边界类型（支持多边界组合，1D暂用单一类型）
        bc_type = cfd.config.boundary_type.lower()
        
        if bc_type == "periodic":
            return PeriodicBoundary(cfd)
        elif bc_type == "dirichlet":
            return DirichletBoundary(cfd)
        elif bc_type == "neumann":
            return NeumannBoundary(cfd)
        else:
            raise ValueError(f"不支持的边界类型：{bc_type}（可选：periodic/dirichlet/neumann）")