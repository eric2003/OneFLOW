# src/infrastructure/__init__.py
"""
基础设施模块
包含网格、域、解、边界条件等基础组件
"""

from .mesh import Mesh
from .domain import Domain
from .solution import Solution
from .boundary import BoundaryCondition, BoundaryConditionFactory

__all__ = [
    'Mesh',
    'Domain', 
    'Solution',
    'BoundaryCondition',
    'BoundaryConditionFactory',
]