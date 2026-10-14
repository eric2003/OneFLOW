param(
    [Parameter(Mandatory = $true)]
    [string]$OneFlowExe,

    [Parameter(Mandatory = $true)]
    [string]$CaseDir,

    [Parameter(Mandatory = $true)]
    [string]$RepositoryRoot,

    [string]$PythonExecutable = "python"
)

$ErrorActionPreference = "Stop"
$resolvedExe = (Resolve-Path -LiteralPath $OneFlowExe).Path
$resolvedCaseDir = (Resolve-Path -LiteralPath $CaseDir).Path
$resolvedRepositoryRoot = (Resolve-Path -LiteralPath $RepositoryRoot).Path
$case = Get-Content -LiteralPath (Join-Path $PSScriptRoot "case.json") -Raw |
    ConvertFrom-Json

$controlFile = Join-Path $resolvedCaseDir "script\control.txt"
$layoutFile = Join-Path $resolvedCaseDir $case.parameters.runtime.layout
$generatedGrid = Join-Path $resolvedCaseDir $case.parameters.inspection.sourceGrid
$oneFlowGrid = Join-Path $resolvedCaseDir $case.parameters.inspection.oneFlowGrid
$gridInfo = Join-Path $resolvedCaseDir $case.parameters.inspection.gridInfo
$packageDir = Join-Path $resolvedCaseDir $case.parameters.inspection.packageDirectory
$exporter = Join-Path $resolvedRepositoryRoot "tools\mesh_export\export_structured_grid.py"

if (-not (Test-Path -LiteralPath $controlFile -PathType Leaf)) {
    throw "Case control file not found: $controlFile"
}
if (-not (Test-Path -LiteralPath $layoutFile -PathType Leaf)) {
    throw "Grid layout file not found: $layoutFile"
}
if (-not (Test-Path -LiteralPath $exporter -PathType Leaf)) {
    throw "Mesh exporter not found: $exporter"
}

$previousLocation = Get-Location
try {
    Set-Location (Split-Path -Parent $resolvedExe)
    & $resolvedExe "r" $resolvedCaseDir
    if ($LASTEXITCODE -ne 0) {
        throw "OneFLOW minimal rectangle case failed with exit code $LASTEXITCODE"
    }
}
finally {
    Set-Location $previousLocation
}

foreach ($outputFile in @($generatedGrid, $oneFlowGrid, $gridInfo)) {
    if (-not (Test-Path -LiteralPath $outputFile -PathType Leaf)) {
        throw "Expected OneFLOW output was not found: $outputFile"
    }
}

& $PythonExecutable $exporter $generatedGrid --output $packageDir --real-bytes 8 --byte-order little
if ($LASTEXITCODE -ne 0) {
    throw "Mesh visualization export failed with exit code $LASTEXITCODE"
}

$summary = Get-Content -LiteralPath (Join-Path $packageDir "summary.json") -Raw |
    ConvertFrom-Json
$expected = $case.parameters.expected
$checks = @(
    [pscustomobject]@{ Name = "zones"; Actual = $summary.counts.zones; Expected = $expected.zones },
    [pscustomobject]@{ Name = "dimension"; Actual = $summary.dimension; Expected = $expected.dimension },
    [pscustomobject]@{ Name = "nodes"; Actual = $summary.counts.nodes; Expected = $expected.nodes },
    [pscustomobject]@{ Name = "cells"; Actual = $summary.counts.cells; Expected = $expected.cells },
    [pscustomobject]@{ Name = "faces"; Actual = $summary.counts.faces; Expected = $expected.faces },
    [pscustomobject]@{ Name = "boundaryFaces"; Actual = $summary.counts.boundaryFaces; Expected = $expected.boundaryFaces }
)
foreach ($check in $checks) {
    if ($check.Actual -ne $check.Expected) {
        throw "Case check failed for $($check.Name): expected $($check.Expected), got $($check.Actual)"
    }
}

$internalFaces = $summary.counts.faces - $summary.counts.boundaryFaces
if ($internalFaces -ne $expected.internalFaces) {
    throw "Case check failed for internalFaces: expected $($expected.internalFaces), got $internalFaces"
}

foreach ($axis in 0..2) {
    if ([Math]::Abs($summary.boundingBox.min[$axis] - $expected.boundingBox.min[$axis]) -gt 1.0e-12) {
        throw "Case check failed for boundingBox.min axis $axis"
    }
    if ([Math]::Abs($summary.boundingBox.max[$axis] - $expected.boundingBox.max[$axis]) -gt 1.0e-12) {
        throw "Case check failed for boundingBox.max axis $axis"
    }
}

Write-Host "Minimal rectangle case passed."
Write-Host "Generated grid: $generatedGrid"
Write-Host "OneFLOW grid: $oneFlowGrid"
Write-Host "Grid summary: $gridInfo"
Write-Host "SVG view: $(Join-Path $packageDir 'views\mesh.svg')"
Write-Host "VTU view: $(Join-Path $packageDir 'views\mesh.vtu')"
