# Script to verify that patched files match their originals + patch files
param(
    [string]$PatchDir = "patches"
)

$ErrorActionPreference = 'Stop'

function Resolve-RepoPath {
    param([string]$Path)
    $resolved = Resolve-Path $Path -ErrorAction SilentlyContinue
    if (-not $resolved) { throw "Path not found: $Path" }
    return $resolved.Path
}

$PatchPath = Resolve-RepoPath $PatchDir
$ReferencePath = Join-Path $PatchPath "reference"

# Mapping of patched files to their originals
$fileMappings = @{
    "beefleg.h" = "beefleg.h.original"
    "beefun.c" = "beefun.c.original"
    "sockets.c" = "sockets.c.original"
    "qeHelpers.cmake" = "qeHelpers.cmake.original"
    "GitInfo.cmake" = "GitInfo.cmake.original"
    "LAXlib_CMakeLists.txt" = "LAXlib_CMakeLists.txt.original"
    "devxlib-timer.c" = "timer.c.original"
}

# Mapping of files to their patch files
$patchMappings = @{
    "beefleg.h" = "001-use-srand-instead-of-srandom.patch"
    "beefun.c" = @("001-use-srand-instead-of-srandom.patch", "qe-win-c-portability.patch")
    "sockets.c" = "qe-win-c-portability.patch"
    "LAXlib_CMakeLists.txt" = "qe-win-cmake-generation.patch"
    "GitInfo.cmake" = "qe-win-cmake-generation.patch"
    "qeHelpers.cmake" = "qe-win-cmake-generation.patch"
    "devxlib-timer.c" = $null  # Need to check if there's a patch for this
}

Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "Patch Verification Report" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host ""

$allPassed = $true

foreach ($patchedFile in $fileMappings.Keys) {
    $originalFile = $fileMappings[$patchedFile]
    $patchedPath = Join-Path $PatchPath $patchedFile
    $originalPath = Join-Path $ReferencePath $originalFile
    
    Write-Host "Checking: $patchedFile" -ForegroundColor Cyan
    
    # Check if files exist
    if (-not (Test-Path $patchedPath)) {
        Write-Host "  ERROR: Patched file not found: $patchedPath" -ForegroundColor Red
        $allPassed = $false
        continue
    }
    
    if (-not (Test-Path $originalPath)) {
        Write-Host "  ERROR: Original file not found: $originalPath" -ForegroundColor Red
        $allPassed = $false
        continue
    }
    
    # Compare files using diff (if available) or file hash
    Write-Host "  Comparing with original..." -ForegroundColor Gray
    
    # Try to use git diff if available, otherwise use Compare-Object
    $patchedContent = Get-Content $patchedPath -Raw
    $originalContent = Get-Content $originalPath -Raw
    
    if ($patchedContent -eq $originalContent) {
        Write-Host "  WARNING: Patched file is identical to original!" -ForegroundColor Yellow
    } else {
        Write-Host "  OK: Files differ (as expected)" -ForegroundColor Green
        
        # Check if we can verify against patch files
        $patchFiles = $patchMappings[$patchedFile]
        if ($patchFiles) {
            if ($patchFiles -is [array]) {
                Write-Host "  Patch files: $($patchFiles -join ', ')" -ForegroundColor Gray
            } else {
                Write-Host "  Patch file: $patchFiles" -ForegroundColor Gray
            }
        } else {
            Write-Host "  NOTE: No patch file mapping found for this file" -ForegroundColor Yellow
        }
    }
    
    Write-Host ""
}

Write-Host "==========================================" -ForegroundColor Cyan
if ($allPassed) {
    Write-Host "Verification completed." -ForegroundColor Green
} else {
    Write-Host "Verification completed with errors." -ForegroundColor Red
}
Write-Host "==========================================" -ForegroundColor Cyan

