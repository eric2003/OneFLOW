@echo off
chcp 65001 >nul

echo ========================================
echo    Step 2: Configuration Physics Update
echo    (Python-based with Intel environment)
echo ========================================
echo.

python --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python not found or not in PATH
    pause
    exit /b 1
)

echo [INFO] Running Python step2 script with full Intel environment support...
echo.

python run_step2.py %*

if errorlevel 1 (
    echo.
    echo [ERROR] Step 2 failed
    pause
    exit /b 1
)

echo.
echo [INFO] Step 2 completed successfully!
echo.
echo [INFO] Next step: Update component manager to support physics
echo [INFO] Run tests independently:
echo       cd build
echo       run_tests.bat
echo.

pause
exit /b 0