# OneFLOW-CFD: 1D Convection Solver

A Fortran implementation of 1D linear convection equation solver with ENO/WENO reconstruction schemes comparison.

## Features

- **Numerical Schemes**: ENO3 and WENO3 reconstruction
- **Time Integration**: Runge-Kutta methods (RK1, RK2)
- **Flux Schemes**: Rusanov and Engquist-Osher fluxes
- **Boundary Conditions**: Periodic boundary conditions
- **Visualization**: Python-based plotting and analysis
- **Automation**: Complete build-run-plot pipeline

## Quick Start

### Prerequisites
- Fortran compiler (Intel ifx/ifort, GNU gfortran, or LLVM flang)
- CMake 3.10+
- Python 3.8+ with NumPy and Matplotlib

### Installation

1. Clone or download the project:
   ```bash
   git clone <repository-url>
   cd OneFLOW-CFD-Fortran