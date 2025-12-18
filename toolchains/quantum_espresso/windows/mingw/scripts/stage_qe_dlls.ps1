# Stage QE DLLs Script
# Copies required DLLs from MinGW to QE installation directory
# Usage: .\scripts\stage_qe_dlls.ps1 [-QeBinDir <path>] [-MingwBinDir <path>]

param(
    [string]$QeBinDir = "",              # QE bin directory (default: auto-detect from script location)
    [string]$MingwBinDir = ""            # MinGW bin directory (default: auto-detect common locations)
)

$ErrorActionPreference = 'Stop'

# Determine script location and repo root
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Resolve-Path (Join-Path $ScriptDir "..\..\..\..\..")

# Set default QE bin directory if not specified
if ([string]::IsNullOrEmpty($QeBinDir)) {
    $QeBinDir = Join-Path $RepoRoot "toolchains\quantum_espresso\windows\mingw\upstream\qe\artifacts\bin"
}

# Find MSYS2/MinGW bin directory if not specified
if ([string]::IsNullOrEmpty($MingwBinDir)) {
    $msys2Paths = @(
        "C:\tools\msys64\mingw64\bin",
        "C:\msys64\mingw64\bin",
        "$env:ProgramFiles\msys64\mingw64\bin"
    )
    
    # Check if MSYS2_ROOT is set
    if ($env:MSYS2_ROOT) {
        $msys2Paths = @(Join-Path $env:MSYS2_ROOT "mingw64\bin") + $msys2Paths
    }
    
    foreach ($path in $msys2Paths) {
        if (Test-Path $path) {
            $MingwBinDir = $path
            break
        }
    }
    
    if ([string]::IsNullOrEmpty($MingwBinDir)) {
        Write-Host "ERROR: MinGW bin directory not found!" -ForegroundColor Red
        Write-Host "Please specify -MingwBinDir parameter or set MSYS2_ROOT environment variable" -ForegroundColor Yellow
        exit 1
    }
}

Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "QE DLL Staging Script" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "QE bin directory: $QeBinDir"
Write-Host "MinGW bin directory: $MingwBinDir"
Write-Host ""

# Verify directories exist
if (-not (Test-Path $QeBinDir)) {
    Write-Host "ERROR: QE bin directory not found: $QeBinDir" -ForegroundColor Red
    exit 1
}

if (-not (Test-Path $MingwBinDir)) {
    Write-Host "ERROR: MinGW bin directory not found: $MingwBinDir" -ForegroundColor Red
    exit 1
}

# Check if pw.exe exists to determine DLL dependencies
$pwExe = Join-Path $QeBinDir "pw.exe"
if (-not (Test-Path $pwExe)) {
    Write-Host "WARNING: pw.exe not found in $QeBinDir" -ForegroundColor Yellow
    Write-Host "Will copy common QE DLLs anyway..." -ForegroundColor Yellow
}

# Required DLLs based on pw.exe dependencies
$requiredDlls = @(
    "libfftw3-3.dll",
    "libfftw3_omp-3.dll",
    "libgcc_s_seh-1.dll",
    "libgfortran-*.dll",          # May have version suffix (e.g., libgfortran-5.dll)
    "libgomp-1.dll",
    "libopenblas*.dll",           # May have version suffix
    "libquadmath-*.dll",          # Required by libgfortran
    "libwinpthread-1.dll"         # Required by libquadmath
)

Write-Host "Copying required DLLs..." -ForegroundColor Yellow
Write-Host ""

$copiedCount = 0
$missingCount = 0
$missingDlls = @()

foreach ($dllPattern in $requiredDlls) {
    $dllFiles = Get-ChildItem -Path $MingwBinDir -Filter $dllPattern -ErrorAction SilentlyContinue
    
    if ($dllFiles) {
        foreach ($dllFile in $dllFiles) {
            $destPath = Join-Path $QeBinDir $dllFile.Name
            try {
                Copy-Item -Path $dllFile.FullName -Destination $destPath -Force
                Write-Host "  Copied: $($dllFile.Name)" -ForegroundColor Green
                $copiedCount++
            } catch {
                Write-Host "  ERROR copying $($dllFile.Name): $_" -ForegroundColor Red
                $missingCount++
            }
        }
    } else {
        Write-Host "  WARNING: Not found: $dllPattern" -ForegroundColor Yellow
        $missingDlls += $dllPattern
        $missingCount++
    }
}

Write-Host ""
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "Summary" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "  DLLs copied: $copiedCount" -ForegroundColor Green
if ($missingCount -gt 0) {
    Write-Host "  Warnings/Missing: $missingCount" -ForegroundColor Yellow
    if ($missingDlls.Count -gt 0) {
        Write-Host "  Missing patterns:" -ForegroundColor Yellow
        foreach ($pattern in $missingDlls) {
            Write-Host "    - $pattern" -ForegroundColor Yellow
        }
    }
}

# Verify critical DLLs are present
Write-Host ""
Write-Host "Verifying critical DLLs..." -ForegroundColor Yellow
$criticalDlls = @("libgomp-1.dll", "libgfortran-*.dll", "libopenblas*.dll", "libfftw3-3.dll")
$allPresent = $true

foreach ($pattern in $criticalDlls) {
    $found = Get-ChildItem -Path $QeBinDir -Filter $pattern -ErrorAction SilentlyContinue
    if ($found) {
        Write-Host "  OK: $pattern" -ForegroundColor Green
    } else {
        Write-Host "  MISSING: $pattern" -ForegroundColor Red
        $allPresent = $false
    }
}

if ($allPresent) {
    Write-Host ""
    Write-Host "All critical DLLs are present!" -ForegroundColor Green
    Write-Host "pw.exe should now run correctly." -ForegroundColor Green
} else {
    Write-Host ""
    Write-Host "WARNING: Some critical DLLs are missing!" -ForegroundColor Yellow
    Write-Host "pw.exe may not run correctly." -ForegroundColor Yellow
}

Write-Host ""

exit 0

