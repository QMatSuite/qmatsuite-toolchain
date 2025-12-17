# QE Source Refresh for Windows oneAPI
param(
    [switch]$DryRun,
    [string]$QeVersion = $null,
    [switch]$NoPatch
)

$ErrorActionPreference = 'Stop'
$ScriptDir   = Split-Path -Parent $MyInvocation.MyCommand.Path
$OneApiRoot  = Split-Path $ScriptDir -Parent   # .../windows/oneapi
$QeUpstreamDir = Join-Path $OneApiRoot "upstream\\qe"
$RepoUrl     = "https://github.com/QEF/q-e.git"

Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "QE Source Refresh (Windows oneAPI)" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "QE target : $QeUpstreamDir"
Write-Host "Repo URL  : $RepoUrl"
if ($DryRun) { Write-Host "DRY RUN - no changes will be made." -ForegroundColor Yellow }
Write-Host ""

function Get-LatestTag {
    $tags = git ls-remote --tags $RepoUrl 2>$null | %{
        if ($_ -match "refs/tags/(qe-[0-9]+\\.[0-9]+(\\.[0-9]+)?)$") { $matches[1] }
    }
    $sorted = $tags | Sort-Object {[Version]($_ -replace "qe-","")} -Descending
    return $sorted[0]
}

if (-not $QeVersion) {
    try {
        $QeVersion = Get-LatestTag
        Write-Host "Selected newest tag: $QeVersion" -ForegroundColor Green
    } catch {
        $QeVersion = "qe-7.5"
        Write-Host "Failed to fetch tags; defaulting to $QeVersion" -ForegroundColor Yellow
    }
} else {
    Write-Host "Using specified tag: $QeVersion" -ForegroundColor Green
}

if (Test-Path $QeUpstreamDir) {
    if ($DryRun) {
        Write-Host "Would remove existing QE directory: $QeUpstreamDir" -ForegroundColor Cyan
    } else {
        Write-Host "Removing existing QE directory: $QeUpstreamDir" -ForegroundColor Yellow
        Remove-Item -Recurse -Force $QeUpstreamDir
    }
}

if (-not $DryRun) {
    New-Item -ItemType Directory -Force -Path (Split-Path $QeUpstreamDir) | Out-Null
    Write-Host "Cloning QE ($QeVersion)..." -ForegroundColor Yellow
    & git clone --recursive --branch $QeVersion --depth 1 $RepoUrl $QeUpstreamDir
    if ($LASTEXITCODE -ne 0) { throw "git clone failed (exit $LASTEXITCODE)" }
    Push-Location $QeUpstreamDir
    try {
        & git submodule update --init --recursive
        if ($LASTEXITCODE -ne 0) { Write-Warning "Submodule init had issues (exit $LASTEXITCODE)" }
    } finally {
        Pop-Location
    }
} else {
    Write-Host "Dry run: would git clone --recursive --branch $QeVersion $RepoUrl $QeUpstreamDir" -ForegroundColor Cyan
}

if (-not $DryRun -and -not $NoPatch) {
    $applyScript = Join-Path $ScriptDir "apply_qe_patches.ps1"
    Write-Host "Applying patches..." -ForegroundColor Yellow
    & powershell -ExecutionPolicy Bypass -File $applyScript -QeDir $QeUpstreamDir
    if ($LASTEXITCODE -ne 0) { throw "Patch application failed (exit $LASTEXITCODE)" }
} elseif ($NoPatch) {
    Write-Host "Skipping patch application (-NoPatch specified)." -ForegroundColor Yellow
}

if (-not $DryRun) {
    $required = @("external/wannier90","external/mbd","external/devxlib")
    $missing = @()
    foreach ($s in $required) {
        if (-not (Test-Path (Join-Path $QeUpstreamDir $s))) { $missing += $s }
    }
    if ($missing.Count -gt 0) {
        Write-Warning "Missing submodules: $($missing -join ', '). Run git submodule update --init --recursive."
    } else {
        Write-Host "Submodules present." -ForegroundColor Green
    }
}

Write-Host "Done. QE clone: $QeUpstreamDir" -ForegroundColor Green

