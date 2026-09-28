# src/numerics/flux/__init__.py
"""
通量计算模块
提供Rusanov、Engquist-Osher等通量计算方法
"""

from .base import InviscidFluxCalculator
from .factory import FluxCalculatorFactory

# 导入具体实现以触发注册
from .rusanov import RusanovFluxCalculator
from .engquist_osher import EngquistOsherFluxCalculator

__all__ = [
    'InviscidFluxCalculator',
    'FluxCalculatorFactory',
    'RusanovFluxCalculator',
    'EngquistOsherFluxCalculator',
]