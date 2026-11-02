# src/boundary/__init__.py
"""
边界条件模块
提供各种边界条件实现和工厂创建接口
"""

from .base import BoundaryCondition
from .factory import BoundaryConditionFactory

# 导入具体类以触发注册
from .periodic import PeriodicBoundary
from .dirichlet import DirichletBoundary
from .neumann import NeumannBoundary

__all__ = [
    'BoundaryCondition',
    'BoundaryConditionFactory',
    'PeriodicBoundary',
    'DirichletBoundary',
    'NeumannBoundary',
]