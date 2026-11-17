# Fortran CFD Project Builder (PowerShell wrapper)

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "  Fortran CFD Project Builder" -ForegroundColor Cyan
Write-Host "  (PowerShell wrapper for Python script)" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""

# 检查Python
$python = Get-Command python -ErrorAction SilentlyContinue
if (-not $python) {
    $python = Get-Command python3 -ErrorAction SilentlyContinue
}

if (-not $python) {
    Write-Host "[ERROR] Python not found. Please install Python 3." -ForegroundColor Red
    Read-Host "Press Enter to exit"
    exit 1
}

# 运行Python构建脚本
Write-Host "[INFO] Running Python build script..." -ForegroundColor Yellow
$argsString = $args -join ' '
$command = "python build.py $argsString"

Invoke-Expression $command

if ($LASTEXITCODE -ne 0) {
    Write-Host "[ERROR] Build failed" -ForegroundColor Red
    Read-Host "Press Enter to exit"
    exit 1
}