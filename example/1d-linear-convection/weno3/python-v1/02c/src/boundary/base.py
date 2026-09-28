# src/boundary/base.py
from abc import ABC, abstractmethod

class BoundaryCondition(ABC):
    """边界条件抽象基类：定义所有边界条件必须实现的接口"""
    def __init__(self, cfd):
        self.cfd = cfd
        self.domain = cfd.domain
        self.config = cfd.config  # 可从配置读取边界参数（如进口速度、固壁温度等）
    
    @abstractmethod
    def apply(self, u):
        """
        应用边界条件到解数组
        :param u: 包含ghost层的解数组（会直接修改该数组）
        :return: None
        """
        pass