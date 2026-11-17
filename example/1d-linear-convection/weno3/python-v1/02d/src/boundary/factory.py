# src/boundary/factory.py
from core.registry import BaseFactory

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
        bc_type = cfd.config.boundary_type
        return BaseFactory.create_component('boundary', bc_type, cfd)
    
    @staticmethod
    def get_available_types():
        """获取所有可用的边界条件类型"""
        return BaseFactory.get_available_components('boundary')