# time_integration/__init__.py

# 导出统一接口
from .factory import TimeIntegratorFactory
from .base import TimeIntegrator

# 触发子模块注册（关键！）
from . import rk1, rk2, rk3