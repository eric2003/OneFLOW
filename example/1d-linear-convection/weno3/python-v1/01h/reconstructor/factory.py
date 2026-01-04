# reconstructor/factory.py
from .eno import EnoReconstructor
#from .weno3 import Weno3Reconstructor
from .weno5 import Weno5Reconstructor
from core.registry import ComponentRegistry, register_component

class ReconstructorFactory:
    @staticmethod
    def create(config, domain):
        scheme = config.recon_scheme.lower()
        order = getattr(config, 'spatial_order', None)
        
        print(f"ReconstructorFactory scheme={scheme}")
        
        if scheme == "weno":
            if order is None:
                raise ValueError("使用 'weno' 时必须设置 config.spatial_order")
            scheme = f"weno{order}"
            
        print(f"ReconstructorFactory scheme={scheme}")
        # 使用BaseFactory，但处理特殊参数
        from core.registry import BaseFactory
        
        print(f"scheme={scheme}")
        if scheme == "eno":
            if order is None:
                order = 3
            # ENO需要特殊参数
            return ComponentRegistry.create('reconstructor', scheme, order, domain.ntcells)
        elif scheme == "weno3" or scheme == "weno5":
            # WENO3,WENO5无参数
            print(f"scheme == 'eno3' or scheme == 'weno5' scheme={scheme}")
            return BaseFactory.create_component('reconstructor', scheme)
        else:
            # 其他情况
            return BaseFactory.create_component('reconstructor', scheme, config)