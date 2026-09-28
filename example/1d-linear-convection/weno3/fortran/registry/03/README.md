# Fortran CFD Project

基于工厂模式和注册系统的CFD求解器框架。

## 项目结构

```
fortran_cfd/
├── CMakeLists.txt              # 根CMake配置
├── scripts/                    # 构建脚本
│   ├── build.py               # Python构建脚本（推荐）
│   ├── build.bat              # Batch包装脚本
│   ├── build.ps1              # PowerShell包装脚本
│   └── requirements.txt       # Python依赖
├── src/                        # 源代码
│   ├── CMakeLists.txt
│   ├── core/                  # 核心模块（注册系统、工厂接口）
│   ├── infrastructure/        # 基础设施（配置、网格）
│   └── numerics/              # 数值方法
├── tests/                      # 测试
│   ├── CMakeLists.txt
│   ├── test_minimal_simple.f90 # 简单测试
│   └── test_factory.f90       # 工厂模式测试
└── README.md                   # 项目说明（本文件）
```

## 快速开始

### 使用Python脚本（推荐）

```bash
# 进入脚本目录
cd scripts

# 基本构建
python build.py

# 清理并构建
python build.py --clean

# 发布构建
python build.py --build-type Release --clean

# 并行构建（使用所有核心）
python build.py -j

# 不运行测试
python build.py --no-tests

# 查看帮助
python build.py --help
```

### 使用批处理脚本

```bash
cd scripts
.\build.bat
```

### 使用PowerShell脚本

```powershell
cd scripts
.\build.ps1
```

## 手动构建

```bash
# 设置Intel环境
cmd.exe "/K" '"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" && powershell'

# 配置项目
mkdir build
cd build
cmake -G "Visual Studio 17 2022" -A x64 -T fortran=ifx ..

# 构建项目
cmake --build . --config Debug

# 运行测试
Debug\test_simple.exe
Debug\test_factory.exe
```

## 功能特性

✅ 组件注册系统  
✅ 工厂模式  
✅ ENO/WENO重构器  
✅ Rusanov通量计算器  
✅ 可扩展架构  
✅ 自动化测试  

## 依赖

- Intel oneAPI（包含ifx编译器）
- CMake 3.12+
- Python 3.6+（仅用于构建脚本）

## 源代码文件

### 核心模块
- `src/core/registry.f90` - 组件注册系统
- `src/core/factory_interfaces.f90` - 工厂接口

### 基础设施
- `src/infrastructure/config.f90` - 配置管理
- `src/infrastructure/mesh.f90` - 网格生成

### 数值方法
- `src/numerics/reconstructor/base.f90` - 重构器基类
- `src/numerics/reconstructor/eno.f90` - ENO重构器
- `src/numerics/reconstructor/weno3.f90` - WENO-3重构器
- `src/numerics/flux/base.f90` - 通量计算器基类
- `src/numerics/flux/rusanov.f90` - Rusanov通量

### 测试
- `tests/test_minimal_simple.f90` - 简单功能测试
- `tests/test_factory.f90` - 工厂模式测试

## 构建脚本

### Python脚本 (build.py)
主构建脚本，支持：
- 自动设置Intel oneAPI环境
- 参数化构建选项
- 并行编译
- 自动化测试
- 彩色输出

### Batch脚本 (build.bat)
Windows批处理包装器，调用Python脚本。

### PowerShell脚本 (build.ps1)
PowerShell包装器，调用Python脚本。

## 示例输出

成功构建后，会看到类似输出：

```
============================================================
  Fortran CFD Project Builder
============================================================

[1/4] Setting up Intel Fortran compiler...
✓ Intel oneAPI environment configured

[2/4] Preparing build directory...
✓ Cleaned build directory

[3/4] Configuring project...
  $ cmake -G Visual Studio 17 2022 -A x64 -T fortran=ifx -DCMAKE_BUILD_TYPE=Debug ..
✓ CMake configuration successful

[4/4] Building project...
  $ cmake --build . --config Debug
✓ Build completed in 12.3 seconds

============================================================
  Running Tests
============================================================

[TEST] Simple functionality test...
✓ Simple test passed

[TEST] Factory pattern test...
✓ Factory test passed

============================================================
  Build Summary
============================================================
✓ Build directory: D:\path\to\build
✓ Build type: Debug
✓ Total time: 15.2s
```

## 开发指南

### 添加新组件

1. 在相应模块中创建Fortran源文件
2. 实现工厂函数
3. 使用`register_component_with_factory`注册组件
4. 添加测试

### 扩展架构

项目采用模块化设计：
1. **注册系统** - 管理所有可用的组件
2. **工厂模式** - 动态创建组件实例
3. **抽象接口** - 确保组件兼容性
4. **配置驱动** - 通过配置文件控制行为

## 故障排除

### 常见问题

1. **Intel环境未找到**
   - 检查Intel oneAPI安装路径
   - 手动运行`setvars.bat`

2. **CMake配置失败**
   - 确保使用正确的生成器：`Visual Studio 17 2022`
   - 检查Fortran编译器是否可用

3. **构建失败**
   - 检查依赖项是否完整
   - 查看详细的错误信息

### 调试建议

- 使用`--verbose`参数获取详细输出
- 检查`build/CMakeCache.txt`了解配置详情
- 查看编译器错误日志

## 许可证

MIT License