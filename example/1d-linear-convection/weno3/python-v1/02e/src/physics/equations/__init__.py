# src/physics/equations/__init__.py
"""
物理方程模块
提供各种控制方程的实现
"""

from .base import Equation
from .linear_advection import LinearAdvectionEquation
from .euler import EulerEquation

__all__ = [
    'Equation',
    'LinearAdvectionEquation',
    'EulerEquation',
]