# src/numerics/time_integration/factory.py
from core.registry import BaseFactory

class TimeIntegratorFactory:
    @staticmethod
    def create(cfd) -> 'TimeIntegrator':
        rk_order = cfd.config.rk_order
        integrator_name = f'rk{rk_order}'
        return BaseFactory.create_component('integrator', integrator_name, cfd)
        
    @staticmethod
    def get_available_types():
        """获取所有可用的重构器类型"""
        from core.registry import BaseFactory
        return BaseFactory.get_available_components('integrator')