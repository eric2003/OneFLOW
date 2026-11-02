param(
    [Parameter(Mandatory = $true)]
    [string]$OneFlowExe,

    [Parameter(Mandatory = $true)]
    [string]$CaseDir
)

$resolvedExe = (Resolve-Path -LiteralPath $OneFlowExe).Path
$resolvedCaseDir = (Resolve-Path -LiteralPath $CaseDir).Path
$caseDir = $resolvedCaseDir
$outputFile = Join-Path $caseDir "results\rae2822_2d.ofl"

$previousLocation = Get-Location
$exitCode = 1
try {
    # Keep runtime resources relative to the executable installation.
    Set-Location (Split-Path -Parent $resolvedExe)
    & $resolvedExe "r" $caseDir
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
Write-Host "OneFLOW grid output: $outputFile"
