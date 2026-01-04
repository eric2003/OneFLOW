# equations/base.py

from abc import ABC, abstractmethod
import numpy as np

class Equation(ABC):
    """
    控制方程抽象基类
    所有物理方程（线性对流、Euler、MHD）必须继承此类
    """

    @property
    @abstractmethod
    def num_equations(self) -> int:
        """返回方程组变量数（标量=1，Euler=3，MHD=8）"""
        pass

    @abstractmethod
    def flux(self, u):
        """
        计算通量函数 f(u)
        :param u: 状态向量 (num_equations,)
        :return: 通量向量 (num_equations,)
        """
        pass

    @abstractmethod
    def max_wave_speed(self, u):
        """
        计算最大波速（用于 CFL 条件）
        :param u: 状态向量
        :return: 标量波速
        """
        pass

    def exact_solution(self, x, t, initial_condition_func):
        """
        可选：提供解析解
        子类可重写；若无解析解，调用方应捕获 NotImplementedError
        :param x: 空间坐标数组 (ncells,)
        :param t: 时间
        :param initial_condition_func: 初值函数 u0(x) -> (num_equations, ncells)
        :return: 解 u(x,t) -> (num_equations, ncells)
        """
        raise NotImplementedError(
            f"Equation '{self.__class__.__name__}' does not provide an exact solution."
        )