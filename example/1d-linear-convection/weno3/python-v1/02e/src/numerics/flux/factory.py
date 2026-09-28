# src/numerics/flux/factory.py
from core.registry import BaseFactory

class FluxCalculatorFactory:
    """通量计算器工厂：隐藏注册类别和组件名称等细节"""

    @staticmethod
    def create(cfd) -> 'InviscidFluxCalculator':
        """
        创建通量计算器实例
        
        Args:
            cfd: CFD 主对象，包含 config.flux_type 等配置
            
        Returns:
            实现 InviscidFluxCalculator 接口的通量计算器实例
        """
        flux_type = cfd.config.flux_type
        return BaseFactory.create_component('flux', flux_type, cfd)
        
    @staticmethod
    def get_available_types():
        """获取所有可用的重构器类型"""
        from core.registry import BaseFactory
        return BaseFactory.get_available_components('flux')        
   