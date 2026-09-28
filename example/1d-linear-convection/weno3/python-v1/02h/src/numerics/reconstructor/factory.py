# src/numerics/reconstructor/factory.py
from core.registry import BaseFactory

class ReconstructorFactory:
    @staticmethod
    def create(cfd):
        config = cfd.config
        scheme = config.recon_scheme.lower()
        
        if scheme.startswith("weno"):
            if scheme == "weno":
                order = getattr(config, 'spatial_order', None)
                scheme = f"weno{order}"
        return BaseFactory.create_component('reconstructor', scheme, cfd)

    @staticmethod
    def get_available_types():
        """获取所有可用的重构器类型"""
        from core.registry import BaseFactory
        return BaseFactory.get_available_components('reconstructor')        