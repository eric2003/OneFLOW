# src/initial_conditions/__init__.py
"""
初始条件模块
提供各种初始条件实现和工厂创建接口
"""

from .base import InitialCondition
from .factory import InitialConditionFactory

from .step import StepFunctionIC
from .sine import SineWaveIC
from .gaussian import GaussianPulseIC
