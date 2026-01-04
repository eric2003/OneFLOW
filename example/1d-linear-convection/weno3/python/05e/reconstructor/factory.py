# reconstructor/factory.py
from .eno import EnoReconstructor
from .weno3 import Weno3Reconstructor

# 导入注册系统
try:
    from cfd_core.registry import ComponentRegistry
except ImportError:
    import sys
    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    try:
        from registry import ComponentRegistry
        print("[reconstructor.factory] 使用本地 registry.py 中的注册系统")
    except ImportError:
        print("[reconstructor.factory] 警告: 未找到注册系统")
        # 提供一个回退方案
        class ComponentRegistry:
            @staticmethod
            def create(category, name, *args, **kwargs):
                if category == 'reconstructor':
                    if name == 'eno':
                        return EnoReconstructor(*args, **kwargs)
                    elif name == 'weno3':
                        return Weno3Reconstructor(*args, **kwargs)
                raise ValueError(f"不支持的 {category}.{name}")

class ReconstructorFactory:
    @staticmethod
    def create(config, domain):
        scheme = config.recon_scheme.lower()
        order = getattr(config, 'spatial_order', None)
        
        if scheme == "weno":
            if order is None:
                raise ValueError("使用 'weno' 时必须设置 config.spatial_order")
            scheme = f"weno{order}"
        
        print(f"[ReconstructorFactory] 请求创建重建器: {scheme}")
        
        try:
            # 根据具体类型传递参数
            if scheme == "eno":
                if order is None:
                    order = 3
                # ENO 需要 order 和 ntcells
                return ComponentRegistry.create('reconstructor', scheme, order, domain.ntcells)
            elif scheme == "weno3":
                # ✅ 关键修复：WENO3 不需要参数！
                return ComponentRegistry.create('reconstructor', scheme)  # 没有参数！
            else:
                # 其他情况默认传递 config
                return ComponentRegistry.create('reconstructor', scheme, config)
                
        except Exception as e:
            print(f"注册系统失败，使用硬编码: {e}")
            # 硬编码回退逻辑
            if scheme == "eno":
                if order is None:
                    order = 3
                return EnoReconstructor(order, domain.ntcells)
            elif scheme == "weno3":
                return Weno3Reconstructor()  # ✅ 无参数！
            else:
                raise ValueError(f"不支持的重建格式：{scheme}")