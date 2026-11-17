# reconstructor/base.py
from abc import ABC, abstractmethod

class Reconstructor(ABC):
    @abstractmethod
    def reconstruct(self, q, cfd):
        pass        