# reconstructor/factory.py
from .eno import EnoReconstructor
from .weno3 import Weno3Reconstructor
from cfd_registry import ComponentRegistry

class ReconstructorFactory:
    @staticmethod
    def create(config, domain):
        scheme = config.recon_scheme.lower()
        order = getattr(config, 'spatial_order', None)
        
        if scheme == "weno":
            if order is None:
                raise ValueError("使用 'weno' 时必须设置 config.spatial_order")
            scheme = f"weno{order}"
        
        # 使用BaseFactory，但处理特殊参数
        from factories.base_factory import BaseFactory
        
        try:
            if scheme == "eno":
                if order is None:
                    order = 3
                # ENO需要特殊参数
                return ComponentRegistry.create('reconstructor', scheme, order, domain.ntcells)
            elif scheme == "weno3":
                # WENO3无参数
                return BaseFactory.create_component('reconstructor', scheme)
            else:
                # 其他情况
                return BaseFactory.create_component('reconstructor', scheme, config)
        except ValueError as e:
            # 简单的回退逻辑
            if scheme == "eno":
                if order is None:
                    order = 3
                from .eno import EnoReconstructor
                return EnoReconstructor(order, domain.ntcells)
            elif scheme == "weno3":
                from .weno3 import Weno3Reconstructor
                return Weno3Reconstructor()
            raise