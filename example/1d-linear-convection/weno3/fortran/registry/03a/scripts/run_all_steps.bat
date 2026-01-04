@echo off
chcp 65001 >nul

echo ========================================
echo    CFD Project: All Steps
echo ========================================
echo.

echo [INFO] Starting Step 1: Physics Modules Test...
call run_step1.bat

if errorlevel 1 (
    echo [ERROR] Step 1 failed
    pause
    exit /b 1
)

echo.
echo [INFO] Starting Step 2: Configuration Physics Update...
call run_step2.bat

if errorlevel 1 (
    echo [ERROR] Step 2 failed
    pause
    exit /b 1
)

echo.
echo ========================================
echo    All Steps Completed Successfully!
echo ========================================
echo.
echo [INFO] Next: Update component manager for physics support
echo [INFO] Run: run_step3.bat (to be created)
echo.

pause
exit /b 0