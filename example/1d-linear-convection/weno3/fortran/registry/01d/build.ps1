# Fortran CFD Project Builder (PowerShell)
Write-Host "=========================================" -ForegroundColor Cyan
Write-Host "  Fortran CFD Project Builder (PowerShell)" -ForegroundColor Cyan
Write-Host "=========================================" -ForegroundColor Cyan
Write-Host ""

# 设置 Intel oneAPI 环境
Write-Host "[1/5] Setting up Intel oneAPI environment..." -ForegroundColor Yellow
$setvarsPath = "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
if (-not (Test-Path $setvarsPath)) {
    Write-Host "[ERROR] Intel oneAPI setvars.bat not found at: $setvarsPath" -ForegroundColor Red
    Write-Host "Please check your Intel oneAPI installation." -ForegroundColor Red
    Read-Host "Press Enter to exit"
    exit 1
}

# 调用 setvars.bat 并设置环境变量
try {
    cmd.exe /c "`"$setvarsPath`" > nul 2>&1 && set" | ForEach-Object {
        if ($_ -match '^([^=]+)=(.*)$') {
            $name = $matches[1]
            $value = $matches[2]
            [Environment]::SetEnvironmentVariable($name, $value, 'Process')
        }
    }
    Write-Host "✓ Intel oneAPI environment set" -ForegroundColor Green
}
catch {
    Write-Host "[ERROR] Failed to set Intel oneAPI environment: $_" -ForegroundColor Red
    Read-Host "Press Enter to exit"
    exit 1
}

# 清理构建目录
Write-Host "[2/5] Cleaning build directory..." -ForegroundColor Yellow
if (Test-Path "build") {
    Remove-Item -Path "build" -Recurse -Force
    Write-Host "✓ Old build directory removed" -ForegroundColor Green
}
else {
    Write-Host "- No existing build directory" -ForegroundColor Gray
}

New-Item -ItemType Directory -Path "build" -Force | Out-Null
Set-Location "build"

# 配置项目
Write-Host "[3/5] Configuring project with Intel Fortran..." -ForegroundColor Yellow
$cmakeArgs = @('-G', 'Visual Studio 17 2022', '-A', 'x64', '-T', 'fortran=ifx', '..')
try {
    & cmake @cmakeArgs
    if ($LASTEXITCODE -ne 0) {
        throw "CMake configuration failed with exit code: $LASTEXITCODE"
    }
    Write-Host "✓ CMake configuration successful" -ForegroundColor Green
}
catch {
    Write-Host "[ERROR] $_" -ForegroundColor Red
    Read-Host "Press Enter to exit"
    exit 1
}

# 构建项目
Write-Host "[4/5] Building project..." -ForegroundColor Yellow
try {
    & cmake --build . --config Debug
    if ($LASTEXITCODE -ne 0) {
        throw "Build failed with exit code: $LASTEXITCODE"
    }
    Write-Host "✓ Build successful" -ForegroundColor Green
}
catch {
    Write-Host "[ERROR] $_" -ForegroundColor Red
    Read-Host "Press Enter to exit"
    exit 1
}

# 运行测试
Write-Host "[5/5] Running tests..." -ForegroundColor Yellow
Write-Host "=========================================" -ForegroundColor Cyan

$testResults = @()

# 运行简单测试
if (Test-Path "Debug\test_simple.exe") {
    Write-Host "[TEST] Simple functionality test..." -ForegroundColor Magenta
    Write-Host "-----------------------------------------" -ForegroundColor DarkGray
    try {
        $output = & "Debug\test_simple.exe"
        Write-Host $output
        $testResults += @{ Name = "Simple Test"; Status = "PASSED" }
        Write-Host ""
    }
    catch {
        Write-Host "[ERROR] Simple test failed: $_" -ForegroundColor Red
        $testResults += @{ Name = "Simple Test"; Status = "FAILED" }
    }
}
else {
    Write-Host "[WARN] test_simple.exe not found" -ForegroundColor Yellow
}

# 运行工厂测试
if (Test-Path "Debug\test_factory.exe") {
    Write-Host "[TEST] Factory pattern test..." -ForegroundColor Magenta
    Write-Host "-----------------------------------------" -ForegroundColor DarkGray
    try {
        $output = & "Debug\test_factory.exe"
        Write-Host $output
        $testResults += @{ Name = "Factory Test"; Status = "PASSED" }
        Write-Host ""
    }
    catch {
        Write-Host "[ERROR] Factory test failed: $_" -ForegroundColor Red
        $testResults += @{ Name = "Factory Test"; Status = "FAILED" }
    }
}
else {
    Write-Host "[WARN] test_factory.exe not found" -ForegroundColor Yellow
}

# 显示测试总结
Write-Host "=========================================" -ForegroundColor Cyan
Write-Host "TEST SUMMARY:" -ForegroundColor Cyan
Write-Host "-----------------------------------------" -ForegroundColor DarkGray
foreach ($test in $testResults) {
    $color = if ($test.Status -eq "PASSED") { "Green" } else { "Red" }
    Write-Host "$($test.Name): $($test.Status)" -ForegroundColor $color
}

Write-Host "=========================================" -ForegroundColor Cyan
Write-Host "Build directory: $(Get-Location)" -ForegroundColor Cyan
Write-Host "=========================================" -ForegroundColor Cyan

if ($testResults.Count -eq 0 -or ($testResults.Status -contains "FAILED")) {
    Read-Host "Press Enter to exit"
}