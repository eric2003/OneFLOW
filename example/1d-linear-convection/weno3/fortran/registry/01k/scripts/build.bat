@echo off
echo ========================================
echo    Fortran CFD Project Builder
echo    (Wrapper for Python script)
echo ========================================
echo.

REM 检查Python
where python >nul 2>nul
if %errorlevel% neq 0 (
    echo [ERROR] Python not found. Please install Python 3.
    pause
    exit /b 1
)

REM 运行Python构建脚本
echo [INFO] Running Python build script...
python build.py %*

if %errorlevel% neq 0 (
    echo [ERROR] Build failed
    pause
    exit /b 1
)

pause