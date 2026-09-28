from abc import ABC, abstractmethod
from core.registry import ComponentRegistry, register_component

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

# ---------------------- 具体边界条件实现（使用装饰器注册）--------------------

@register_component('boundary', 'periodic')
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
            
@register_component('boundary', 'dirichlet')
class DirichletBoundary(BoundaryCondition):
    """Dirichlet（固定值）边界条件（如进口固定速度、固壁零速度）"""
    def apply(self, u):
        nghosts = self.domain.nghosts
        ist = self.domain.ist
        ied = self.domain.ied
        
        # ✅ 修复：使用 getattr 而不是 .get() 方法
        # 旧代码: left_value = self.config.get("left_boundary_value", 1.0)
        # 新代码: 使用 getattr，它对于类和实例都适用
        left_value = getattr(self.config, "left_boundary_value", 1.0)
        
        # 左边界（进口）固定值
        for ig in range(nghosts):
            u[ist - 1 - ig] = left_value
        
        # ✅ 修复：同样使用 getattr
        right_value = getattr(self.config, "right_boundary_value", 2.0)
        
        # 右边界（出口）固定值
        for ig in range(nghosts):
            u[ied + ig] = right_value
        
        # 调试信息
        if hasattr(self.config, 'debug') and self.config.debug:
            print(f"  应用Dirichlet边界: 左值={left_value}, 右值={right_value}")
            
@register_component('boundary', 'neumann')
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
        
        # 调试信息
        if hasattr(self.config, 'debug') and self.config.debug:
            print(f"  应用Neumann边界: 零梯度")

class BoundaryConditionFactory:
    """边界条件工厂"""
    
    @staticmethod
    def create(cfd) -> 'BoundaryCondition':
        """
        创建边界条件实例
        
        Args:
            cfd: CFD对象
        
        Returns:
            边界条件实例
        """
        from core.registry import BaseFactory
        bc_type = cfd.config.boundary_type
        return BaseFactory.create_component('boundary', bc_type, cfd)
        
