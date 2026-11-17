# src/physics/systems/__init__.py
"""
物理系统模块
组合方程和解析解逻辑
"""

from .base import PhysicalSystem
from .factory import PhysicalSystemFactory
from .linear_advection_system import LinearAdvectionSystem

__all__ = [
    'PhysicalSystem',
    'PhysicalSystemFactory',
    'LinearAdvectionSystem',
]