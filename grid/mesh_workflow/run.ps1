param(
    [Parameter(Mandatory = $true)]
    [string]$OneFlowExe
)

$resolvedExe = (Resolve-Path -LiteralPath $OneFlowExe).Path
$caseDir = $PSScriptRoot
$outputDir = Join-Path $caseDir "results"
New-Item -ItemType Directory -Force -Path $outputDir | Out-Null

$previousLocation = Get-Location
$exitCode = 1
try {
    Set-Location (Split-Path -Parent $resolvedExe)
    & $resolvedExe 0 $caseDir
    $exitCode = $LASTEXITCODE
}
finally {
    Set-Location $previousLocation
}

if ($exitCode -ne 0) {
    Write-Error "OneFLOW grid workflow failed with exit code $exitCode"
    exit $exitCode
}

Write-Host "Mesh workflow completed."
Write-Host "OneFLOW grid output: $(Join-Path $outputDir 'rae2822_2d.ofl')"
