# src/physics/problems/__init__.py
"""
问题定义模块
"""

from .base import Problem
from .factory import ProblemFactory
from .linear_advection import LinearAdvectionProblem

__all__ = [
    'Problem',
    'ProblemFactory',
    'LinearAdvectionProblem',
]