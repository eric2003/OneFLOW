@echo off
echo =========================================
echo     Fortran CFD Project Builder (Batch)
echo =========================================
echo.

REM 设置 Intel oneAPI 环境
echo [1/5] Setting up Intel oneAPI environment...
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" > nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Failed to set Intel oneAPI environment
    pause
    exit /b 1
)
echo ✓ Intel oneAPI environment set

REM 清理构建目录
echo [2/5] Cleaning build directory...
if exist build (
    rmdir /s /q build
    echo ✓ Old build directory removed
) else (
    echo - No existing build directory
)
mkdir build
cd build

REM 配置项目
echo [3/5] Configuring project with Intel Fortran...
cmake -G "Visual Studio 17 2022" -A x64 -T fortran=ifx ..
if %errorlevel% neq 0 (
    echo [ERROR] CMake configuration failed
    pause
    exit /b 1
)
echo ✓ CMake configuration successful

REM 构建项目
echo [4/5] Building project...
cmake --build . --config Debug
if %errorlevel% neq 0 (
    echo [ERROR] Build failed
    pause
    exit /b 1
)
echo ✓ Build successful

REM 运行测试
echo [5/5] Running tests...
echo =========================================
if exist "Debug\test_simple.exe" (
    echo [TEST] Simple functionality test...
    echo -----------------------------------------
    Debug\test_simple.exe
    echo.
) else (
    echo [WARN] test_simple.exe not found
)

if exist "Debug\test_factory.exe" (
    echo [TEST] Factory pattern test...
    echo -----------------------------------------
    Debug\test_factory.exe
    echo.
) else (
    echo [WARN] test_factory.exe not found
)

echo =========================================
echo Build directory: %CD%
echo =========================================
pause