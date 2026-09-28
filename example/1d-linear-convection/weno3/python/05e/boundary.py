from abc import ABC, abstractmethod

# ---------------------- 导入注册系统 ----------------------
try:
    # 首先尝试从新的核心位置导入
    from cfd_core.registry import ComponentRegistry, register_component
except ImportError:
    # 如果新的核心模块还不存在，使用我们之前定义的简化版本
    # 这将确保代码立即可用，未来再统一迁移到核心模块
    import sys
    import os
    # 添加当前目录，以便找到可能在同一文件夹的 registry.py
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    try:
        from registry import ComponentRegistry, register_component
        print("[boundary] 使用本地 registry.py 中的注册系统")
    except ImportError:
        # 如果完全找不到，提供一个极简定义防止报错（仅用于过渡）
        print("[boundary] 警告: 未找到注册系统，使用极简回退方案")
        class ComponentRegistry:
            _registries = {}
            @classmethod
            def register(cls, *args, **kwargs): pass
            @classmethod
            def create(cls, *args, **kwargs):
                # 这是一个非常基本的回退，仅用于演示
                # 在实际替换前，您应确保 registry.py 存在
                raise RuntimeError("注册系统未正确初始化。请确保 registry.py 在 Python 路径中。")
        def register_component(*args, **kwargs):
            def dummy_decorator(cls):
                return cls
            return dummy_decorator

# ---------------------- 边界条件抽象基类（统一接口） ----------------------
class BoundaryCondition(ABC):
    """边界条件抽象基类：定义所有边界条件必须实现的接口"""
    def __init__(self, cfd):
        self.cfd = cfd
        self.domain = cfd.domain
        self.config = cfd.config  # 可从配置读取边界参数（如进口速度、固壁温度等）
    
    @abstractmethod
    def apply(self, u):
        """
        应用边界条件到解数组
        :param u: 包含ghost层的解数组（会直接修改该数组）
        :return: None
        """
        pass

# ---------------------- 具体边界条件实现（使用装饰器注册）--------------------

@register_component('boundary', 'periodic')
class PeriodicBoundary(BoundaryCondition):
    """周期边界条件（1D专用）"""
    def apply(self, u):
        nghosts = self.domain.nghosts
        ist = self.domain.ist
        ied = self.domain.ied
        
        # 左ghost层 = 右物理层
        for ig in range(nghosts):
            u[ist - 1 - ig] = u[ied - 1 - ig]
        
        # 右ghost层 = 左物理层
        for ig in range(nghosts):
            u[ied + ig] = u[ist + ig]
            
@register_component('boundary', 'dirichlet')
class DirichletBoundary(BoundaryCondition):
    """Dirichlet（固定值）边界条件（如进口固定速度、固壁零速度）"""
    def apply(self, u):
        nghosts = self.domain.nghosts
        ist = self.domain.ist
        ied = self.domain.ied
        
        # 左边界（进口）固定值（从配置读取）
        left_value = self.config.get("left_boundary_value", 1.0)
        for ig in range(nghosts):
            u[ist - 1 - ig] = left_value
        
        # 右边界（出口）固定值（从配置读取）
        right_value = self.config.get("right_boundary_value", 2.0)
        for ig in range(nghosts):
            u[ied + ig] = right_value
            
@register_component('boundary', 'neumann')
class NeumannBoundary(BoundaryCondition):
    """Neumann（零梯度）边界条件（如出口无梯度）"""
    def apply(self, u):
        nghosts = self.domain.nghosts
        ist = self.domain.ist
        ied = self.domain.ied
        
        # 左边界零梯度
        for ig in range(nghosts):
            u[ist - 1 - ig] = u[ist + ig]
        
        # 右边界零梯度
        for ig in range(nghosts):
            u[ied + ig] = u[ied - 1 - ig]

# ---------------------- 新的边界条件工厂（使用注册系统） ----------------------
class BoundaryConditionFactory:
    """边界条件工厂：根据配置创建对应边界条件实例（新版本，使用注册系统）"""
    
    @staticmethod
    def create(cfd):
        """
        使用注册系统创建边界条件实例。
        保持与旧版本完全相同的接口，实现无缝替换。
        """
        # 从配置读取边界类型
        bc_type = cfd.config.boundary_type.lower()
        
        try:
            # 使用注册系统创建实例
            # ComponentRegistry.create(类别, 名称, 传递给构造函数的参数)
            bc_instance = ComponentRegistry.create('boundary', bc_type, cfd)
            return bc_instance
        except (ValueError, RuntimeError) as e:
            # 增强错误信息，列出可用的选项
            available = []
            try:
                # 尝试从注册表获取可用列表
                all_comps = ComponentRegistry.list_all()
                available = all_comps.get('boundary', [])
            except:
                # 如果注册表不可用，使用硬编码列表作为后备
                available = ['periodic', 'dirichlet', 'neumann']
            
            raise ValueError(
                f"不支持的边界类型：'{bc_type}'\n"
                f"可用类型：{available}\n"
                f"原始错误：{e}"
            )