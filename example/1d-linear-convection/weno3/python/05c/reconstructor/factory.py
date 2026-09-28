# reconstructor/factory.py
from .eno import EnoReconstructor
from .weno3 import Weno3Reconstructor

class ReconstructorFactory:
    _schemes = {
        "eno": EnoReconstructor,
        "weno3": Weno3Reconstructor,      # ← 注意：这里用 "weno3"
        # "weno5": Weno5Reconstructor,   # 未来扩展
    }

    @staticmethod
    def create(config, domain):
        scheme = config.recon_scheme.lower()
        order = getattr(config, 'spatial_order', None)
        
        # ✅ 关键：将 "weno" 自动映射为 "weno3", "weno5" 等
        if scheme == "weno":
            if order is None:
                raise ValueError("使用 'weno' 时必须设置 config.spatial_order")
            scheme = f"weno{order}"

        # 检查是否支持
        if scheme not in ReconstructorFactory._schemes:
            supported = list(ReconstructorFactory._schemes.keys())
            raise ValueError(f"不支持的重建格式：'{scheme}'（支持：{supported}）")

        recon_cls = ReconstructorFactory._schemes[scheme]

        # 根据 scheme 类型创建实例
        if scheme.startswith("eno"):
            return recon_cls(order, domain.ntcells)
        elif scheme.startswith("weno"):
            return recon_cls()  # WENO 类通常无参
        # 可扩展 elif...
