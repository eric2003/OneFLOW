# reconstructor/base.py
from abc import ABC, abstractmethod

class Reconstructor(ABC):
    @abstractmethod
    def compute_face_values(self, q, cfd):
        pass        