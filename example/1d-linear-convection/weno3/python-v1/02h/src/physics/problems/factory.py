# src/physics/problems/factory.py
"""
问题工厂
"""

from core.registry import BaseFactory

class ProblemFactory:
    """问题工厂"""
    
    @staticmethod
    def create(problem_type: str, config):
        """
        创建问题实例
        
        Args:
            problem_type: 问题类型，如 'linear_advection'
            config: 配置对象
        
        Returns:
            问题实例
        """
        return BaseFactory.create_component('problem', problem_type, config)
    
    @staticmethod
    def get_available_types():
        """获取所有可用的问题类型"""
        return BaseFactory.get_available_components('problem')