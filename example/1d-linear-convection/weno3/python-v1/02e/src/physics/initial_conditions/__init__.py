# src/physics/initial_conditions/__init__.py
from .base import InitialCondition
from .factory import InitialConditionFactory

# 导入具体类以触发注册
from .step import StepFunctionIC
from .sine import SineWaveIC
from .gaussian import GaussianPulseIC

__all__ = [
    'InitialCondition',
    'InitialConditionFactory',
    'StepFunctionIC',
    'SineWaveIC',
    'GaussianPulseIC',
]