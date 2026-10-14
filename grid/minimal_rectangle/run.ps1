param(
    [Parameter(Mandatory = $true)]
    [string]$TestExecutable,

    [Parameter(Mandatory = $true)]
    [string]$WorkspaceDir,

    [string]$PythonExecutable = "python"
)

$ErrorActionPreference = "Stop"
$case = Get-Content -LiteralPath (Join-Path $PSScriptRoot "case.json") -Raw |
    ConvertFrom-Json
$repositoryRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$testExecutablePath = (Resolve-Path -LiteralPath $TestExecutable).Path
$workspacePath = [System.IO.Path]::GetFullPath($WorkspaceDir)
$generatedDir = Join-Path $workspacePath "generated"
$packageDir = Join-Path $workspacePath $case.parameters.inspection.packageDirectory
$gridPath = Join-Path $generatedDir "minimalRectangle.ofl"
$exporterPath = Join-Path $repositoryRoot "tools\mesh_export\export_structured_grid.py"

if ($case.parameters.generation.dimension -ne 2) {
    throw "This case currently supports dimension 2 only."
}

New-Item -ItemType Directory -Force -Path $workspacePath | Out-Null

$previousOutputDir = $env:ONEFLOW_MESH_CASE_OUTPUT_DIR
$previousLayoutFile = $env:ONEFLOW_MESH_CASE_LAYOUT_FILE
$previousKeepOutput = $env:ONEFLOW_KEEP_MESH_EXPORT_CASE
try {
    $env:ONEFLOW_MESH_CASE_OUTPUT_DIR = $generatedDir
    $env:ONEFLOW_MESH_CASE_LAYOUT_FILE = Join-Path $repositoryRoot $case.parameters.generation.layout
    $env:ONEFLOW_KEEP_MESH_EXPORT_CASE = "1"

    $filter = "--gtest_filter=$($case.parameters.generation.testFilter)"
    & $testExecutablePath $filter
    if ($LASTEXITCODE -ne 0) {
        throw "Minimal rectangle generation test failed with exit code $LASTEXITCODE"
    }
    if (-not (Test-Path -LiteralPath $gridPath -PathType Leaf)) {
        throw "Generated grid file was not found: $gridPath"
    }

    & $PythonExecutable $exporterPath $gridPath `
        --output $packageDir `
        --real-bytes $case.parameters.inspection.realBytes `
        --byte-order $case.parameters.inspection.byteOrder
    if ($LASTEXITCODE -ne 0) {
        throw "Mesh export failed with exit code $LASTEXITCODE"
    }

    $summary = Get-Content -LiteralPath (Join-Path $packageDir "summary.json") -Raw |
        ConvertFrom-Json
    $expected = $case.parameters.expected
    $checks = @(
        [pscustomobject]@{ Name = "zones"; Actual = $summary.counts.zones; Expected = $expected.zones },
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
    Write-Host "Generated grid: $gridPath"
    Write-Host "Mesh package: $packageDir"
    Write-Host "SVG view: $(Join-Path $packageDir 'views\mesh.svg')"
    Write-Host "VTU view: $(Join-Path $packageDir 'views\mesh.vtu')"
}
finally {
    $env:ONEFLOW_MESH_CASE_OUTPUT_DIR = $previousOutputDir
    $env:ONEFLOW_MESH_CASE_LAYOUT_FILE = $previousLayoutFile
    $env:ONEFLOW_KEEP_MESH_EXPORT_CASE = $previousKeepOutput
}
