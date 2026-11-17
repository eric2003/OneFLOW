# src/numerics/time_integration/__init__.py
"""
时间积分模块
提供显式Runge-Kutta方法
"""

from .base import TimeIntegrator
from .factory import TimeIntegratorFactory

# 导入具体实现以触发注册
from .rk1 import RK1Integrator
from .rk2 import RK2Integrator
from .rk3 import RK3Integrator

__all__ = [
    'TimeIntegrator',
    'TimeIntegratorFactory',
    'RK1Integrator',
    'RK2Integrator',
    'RK3Integrator',
]