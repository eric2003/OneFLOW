@echo off
chcp 65001 >nul

echo ========================================
echo    Step 1: Physics Modules Test
echo    (Python-based with Intel environment)
echo ========================================
echo.

python --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python not found or not in PATH
    pause
    exit /b 1
)

echo [INFO] Running Python step1 script with full Intel environment support...
echo.

python run_step1.py %*

if errorlevel 1 (
    echo.
    echo [ERROR] Step 1 failed
    pause
    exit /b 1
)

echo.
echo [INFO] Step 1 completed successfully!
echo.
echo [INFO] Next step: Update config to include physics settings
echo [INFO] Run tests independently:
echo       cd build
echo       run_tests.bat
echo.

pause
exit /b 0