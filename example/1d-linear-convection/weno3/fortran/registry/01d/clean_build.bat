@echo off
echo === Clean Build for Fortran CFD ===
echo.

REM Clean build directory
if exist build rmdir /s /q build
if exist modules rmdir /s /q modules
mkdir build
cd build

echo 1. Configuring project...
cmake -G "Visual Studio 17 2022" -A x64 -T fortran=ifx ..

if errorlevel 1 (
    echo Configuration failed!
    pause
    exit /b 1
)

echo.
echo 2. Building project...
cmake --build . --config Debug

if errorlevel 1 (
    echo Build failed!
    pause
    exit /b 1
)

echo.
echo 3. Running tests...
if exist "Debug\test_simple.exe" (
    echo "=== Simple Test ==="
    Debug\test_simple.exe
    echo.
)

if exist "Debug\test_factory.exe" (
    echo "=== Factory Test ==="
    Debug\test_factory.exe
) else (
    echo Test executable not found!
)

echo.
echo === Build process completed ===
echo Build directory: %CD%
pause