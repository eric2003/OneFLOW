# mesh.py
import numpy as np

# Mesh class: defines computational grid
class Mesh:
    def __init__(self):
        self.xmin = 0.0
        self.xmax = 2.0
        self.ncells = 40
        self.nnodes = self.ncells + 1
        self.nx = self.ncells
        self.x = np.zeros(self.nnodes)
        self.xcc = np.zeros(self.ncells)
        self.init_mesh()
        
    def init_mesh(self):
        self.L = self.xmax - self.xmin
        self.dx = self.L / self.ncells
        
        # Generate node coordinates
        for i in range(self.nnodes):
            self.x[i] = self.xmin + i * self.dx
            
        # Generate cell center coordinates
        for i in range(self.ncells):
            self.xcc[i] = 0.5 * (self.x[i] + self.x[i+1])