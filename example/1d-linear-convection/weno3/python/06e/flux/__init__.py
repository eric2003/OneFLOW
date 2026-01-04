# flux/__init__.py
from .base import InviscidFluxCalculator
from .factory import FluxCalculatorFactory

# 确保子模块被导入以触发注册
from . import rusanov
from . import engquist_osher