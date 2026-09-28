# OneFLOW-CFD C++: 1D Convection Solver

A C++ implementation of 1D linear convection equation solver with ENO/WENO reconstruction schemes comparison, converted from Fortran version.

## Features

- **Numerical Schemes**: ENO3 and WENO3 reconstruction
- **Time Integration**: Runge-Kutta methods (RK1, RK2)
- **Flux Schemes**: Rusanov and Engquist-Osher fluxes
- **Boundary Conditions**: Periodic boundary conditions
- **Visualization**: Python-based plotting and analysis
- **Automation**: Complete build-run-plot pipeline
- **Modern C++**: C++17 with RAII, smart pointers, and STL containers

## Project Structure

OneFLOW-CFD-CPP/
├── CMakeLists.txt # CMake构建配置
├── include/
│ └── cfd_solver.hpp # 头文件（类声明）
├── src/
│ ├── cfd_solver.cpp # 实现文件
│ └── main.cpp # 主程序
├── scripts/
│ ├── postprocess.py # 后处理Python脚本
│ └── run_all.py # 自动化脚本
└── results/ # 输出结果目录


## Key C++ Features

1. **Modern C++17**: Using `std::vector`, `std::unique_ptr`, RAII
2. **Object-Oriented Design**: Abstract base class for reconstructors
3. **Namespace Organization**: All CFD code in `cfd` namespace
4. **Error Handling**: Using exceptions and runtime errors
5. **Memory Safety**: Smart pointers and RAII for automatic cleanup
6. **Type Safety**: Strongly typed using `using` aliases

## Quick Start

### Prerequisites
- C++ compiler with C++17 support (GCC 7+, Clang 5+, MSVC 2017+)
- CMake 3.10+
- Python 3.8+ with NumPy and Matplotlib

### Build and Run

**Method 1: Manual Build**
```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release
./oneflow_cfd
python ../scripts/postprocess.py
```

**Method 2: Automated Pipeline


```bash
# Make script executable (Linux/Mac)
chmod +x scripts/run_all.py

# Run complete automation
python scripts/run_all.py
```

**Method 3: With Options

```bash
python scripts/run_all.py --verbose  # Detailed output
python scripts/run_all.py --build Debug  # Debug build
python scripts/run_all.py --clean    # Clean before building
```