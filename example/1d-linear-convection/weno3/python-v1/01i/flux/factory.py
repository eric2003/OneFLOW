# flux/factory.py
"""
通量计算器专用工厂（封装注册细节，提供清晰接口）
符合你希望“将创建逻辑封装在工厂中”的设计原则
"""

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
   