@echo off
chcp 65001 >nul

echo ========================================
echo    Fortran CFD Project Builder
echo    (Python-based with Intel environment)
echo ========================================
echo.

python --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python not found or not in PATH
    pause
    exit /b 1
)

echo [INFO] Running Python build script with full Intel environment support...
echo.

python build.py %*

if errorlevel 1 (
    echo.
    echo [ERROR] Build failed
    pause
    exit /b 1
)

echo.
echo [INFO] Build completed successfully!
echo.
echo [INFO] To run tests independently:
echo       cd build
echo       run_tests.bat
echo.

pause
exit /b 0