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
$runner = Join-Path $resolvedRepositoryRoot "tools\mesh_workflow.py"
$caseFile = Join-Path $PSScriptRoot "case.json"

if (-not (Test-Path -LiteralPath $runner -PathType Leaf)) {
    throw "Mesh workflow runner not found: $runner"
}

& $PythonExecutable $runner run `
    --case $caseFile `
    --workdir $resolvedCaseDir `
    --oneflow $resolvedExe `
    --repository-root $resolvedRepositoryRoot

if ($LASTEXITCODE -ne 0) {
    throw "Mesh workflow failed with exit code $LASTEXITCODE"
}
