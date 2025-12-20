# Stage QE Executables and DLLs Script
# Copies executables from artifacts and required DLLs from MSYS2 UCRT64 to distribution directory
# Usage: .\scripts\stage_qe_dlls.ps1 [-ArtifactsDir <path>] [-DistDir <path>] [-MingwBinDir <path>]

param(
    [string]$ArtifactsDir = "",          # QE artifacts directory (default: auto-detect)
    [string]$DistDir = "",                # Distribution directory (default: toolchains\quantum_espresso\windows\mingw\dist\win-mingw)
    [string]$MingwBinDir = ""            # MSYS2 UCRT64 bin directory (default: auto-detect)
)

$ErrorActionPreference = 'Stop'

# Determine script location and base directory (mingw directory)
$ScriptDir = $PSScriptRoot
if (-not $ScriptDir) {
    $ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
}

# Script is in scripts\ subdirectory, so base directory is parent
$MingwBaseDir = Split-Path -Parent $ScriptDir

# Set default artifacts directory if not specified (relative to script location)
if ([string]::IsNullOrEmpty($ArtifactsDir)) {
    $ArtifactsDir = Join-Path $MingwBaseDir "upstream\qe\artifacts"
}

# Set default distribution directory if not specified (relative to script location)
if ([string]::IsNullOrEmpty($DistDir)) {
    $DistDir = Join-Path $MingwBaseDir "dist\win-mingw"
}

# Resolve to absolute paths (relative to script location, not working directory)
# First, ensure MingwBaseDir is absolute
try {
    $MingwBaseDir = (Resolve-Path $MingwBaseDir).Path
} catch {
    # If it doesn't exist, make it absolute based on script location
    $MingwBaseDir = [System.IO.Path]::GetFullPath($MingwBaseDir)
}

# Resolve ArtifactsDir relative to script base directory
if (-not [System.IO.Path]::IsPathRooted($ArtifactsDir)) {
    # If relative, resolve relative to script base directory
    $ArtifactsDir = Join-Path $MingwBaseDir $ArtifactsDir
}

# Resolve DistDir relative to script base directory
if (-not [System.IO.Path]::IsPathRooted($DistDir)) {
    # If relative, resolve relative to script base directory
    $DistDir = Join-Path $MingwBaseDir $DistDir
}

# Convert to fully qualified absolute paths
try {
    $ArtifactsDir = [System.IO.Path]::GetFullPath($ArtifactsDir)
} catch {
    # If GetFullPath fails, construct absolute path manually
    if (-not [System.IO.Path]::IsPathRooted($ArtifactsDir)) {
        $ArtifactsDir = [System.IO.Path]::Combine($MingwBaseDir, $ArtifactsDir)
    }
    $ArtifactsDir = [System.IO.Path]::GetFullPath($ArtifactsDir)
}

try {
    $DistDir = [System.IO.Path]::GetFullPath($DistDir)
} catch {
    # If GetFullPath fails, construct absolute path manually
    if (-not [System.IO.Path]::IsPathRooted($DistDir)) {
        $DistDir = [System.IO.Path]::Combine($MingwBaseDir, $DistDir)
    }
    $DistDir = [System.IO.Path]::GetFullPath($DistDir)
}

# Create dist directory if it doesn't exist
if (-not (Test-Path $DistDir)) {
    try {
        $created = New-Item -ItemType Directory -Path $DistDir -Force -ErrorAction Stop
        $DistDir = $created.FullName
        Write-Host "Created distribution directory: $DistDir" -ForegroundColor Green
    } catch {
        Write-Host "ERROR: Failed to create distribution directory: $DistDir" -ForegroundColor Red
        Write-Host "Error: $_" -ForegroundColor Red
        exit 1
    }
} else {
    try {
        $DistDir = (Resolve-Path $DistDir).Path
    } catch {
        # Empty directory, use path as-is
    }
}

# Find MSYS2 UCRT64 bin directory if not specified
if ([string]::IsNullOrEmpty($MingwBinDir)) {
    $msys2Paths = @(
        "C:\msys64\ucrt64\bin",
        "C:\tools\msys64\ucrt64\bin",
        "$env:ProgramFiles\msys64\ucrt64\bin"
    )
    
    # Check if MSYS2_ROOT is set
    if ($env:MSYS2_ROOT) {
        $msys2Paths = @(Join-Path $env:MSYS2_ROOT "ucrt64\bin") + $msys2Paths
    }
    
    foreach ($path in $msys2Paths) {
        if (Test-Path $path) {
            $MingwBinDir = $path
            break
        }
    }
    
    if ([string]::IsNullOrEmpty($MingwBinDir)) {
        Write-Host "ERROR: MSYS2 UCRT64 bin directory not found!" -ForegroundColor Red
        Write-Host "Please specify -MingwBinDir parameter or set MSYS2_ROOT environment variable" -ForegroundColor Yellow
        Write-Host "Expected locations:" -ForegroundColor Yellow
        foreach ($path in $msys2Paths) {
            Write-Host "  - $path" -ForegroundColor Gray
        }
        exit 1
    }
}

Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "QE Staging Script" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Source (artifacts) directory: $ArtifactsDir" -ForegroundColor White
Write-Host "Destination (distribution) directory: $DistDir" -ForegroundColor White
Write-Host "MSYS2 UCRT64 bin directory: $MingwBinDir" -ForegroundColor White
Write-Host ""

# Verify directories exist
if (-not (Test-Path $ArtifactsDir)) {
    Write-Host "ERROR: Artifacts directory not found: $ArtifactsDir" -ForegroundColor Red
    exit 1
}

if (-not (Test-Path $MingwBinDir)) {
    Write-Host "ERROR: MSYS2 UCRT64 bin directory not found: $MingwBinDir" -ForegroundColor Red
    exit 1
}

# Step 1: Clean distribution directory (remove all files)
Write-Host "Cleaning distribution directory..." -ForegroundColor Yellow
$allFiles = Get-ChildItem -Path $DistDir -File -ErrorAction SilentlyContinue
if ($allFiles) {
    foreach ($file in $allFiles) {
        Remove-Item -Path $file.FullName -Force
    }
    Write-Host "  Removed $($allFiles.Count) existing file(s)" -ForegroundColor Cyan
} else {
    Write-Host "  Distribution directory is empty" -ForegroundColor Gray
}
Write-Host ""

# Step 2: Find all *.x executables recursively in artifacts directory
Write-Host "Finding executables in artifacts directory..." -ForegroundColor Yellow
Write-Host "Searching: $ArtifactsDir" -ForegroundColor Gray
Write-Host ""

$xFiles = Get-ChildItem -Path $ArtifactsDir -Filter "*.x" -Recurse -File -ErrorAction SilentlyContinue

if ($null -eq $xFiles -or $xFiles.Count -eq 0) {
    Write-Host "  WARNING: No .x files found in $ArtifactsDir" -ForegroundColor Yellow
    Write-Host "  Make sure QE has been built and installed." -ForegroundColor Yellow
} else {
    Write-Host "  Found $($xFiles.Count) .x file(s)" -ForegroundColor Cyan
    foreach ($file in $xFiles) {
        $relativePath = $file.FullName.Substring($ArtifactsDir.Length + 1)
        Write-Host "    - $relativePath" -ForegroundColor Gray
    }
}
Write-Host ""

# Step 3: Copy all *.x files to distribution directory (flat layout) and rename to *.exe
Write-Host "Copying executables to distribution directory..." -ForegroundColor Yellow
Write-Host "Destination: $DistDir" -ForegroundColor Cyan
Write-Host ""

$copiedExeCount = 0
$skippedCount = 0

if ($xFiles -and $xFiles.Count -gt 0) {
    foreach ($xFile in $xFiles) {
        $baseName = [System.IO.Path]::GetFileNameWithoutExtension($xFile.Name)
        $exeName = "$baseName.exe"
        $destPath = Join-Path $DistDir $exeName
        
        # Check for name conflicts (multiple *.x files with same base name)
        if (Test-Path $destPath) {
            Write-Host "  WARNING: Skipping $($xFile.Name) - $exeName already exists" -ForegroundColor Yellow
            $skippedCount++
            continue
        }
        
        try {
            # Verify destination directory exists
            if (-not (Test-Path $DistDir)) {
                Write-Host "  ERROR: Destination directory does not exist: $DistDir" -ForegroundColor Red
                exit 1
            }
            Copy-Item -Path $xFile.FullName -Destination $destPath -Force
            Write-Host "  Copied: $($xFile.Name) -> $exeName" -ForegroundColor Green
            Write-Host "    To: $destPath" -ForegroundColor Gray
            $copiedExeCount++
        } catch {
            Write-Host "  ERROR copying $($xFile.Name): $_" -ForegroundColor Red
            Write-Host "    Source: $($xFile.FullName)" -ForegroundColor Gray
            Write-Host "    Destination: $destPath" -ForegroundColor Gray
        }
    }
    
    Write-Host ""
    Write-Host "Copied $copiedExeCount executable(s)" -ForegroundColor Cyan
    if ($skippedCount -gt 0) {
        Write-Host "Skipped $skippedCount executable(s) due to name conflicts" -ForegroundColor Yellow
    }
    Write-Host ""
} else {
    Write-Host "  No executables to copy" -ForegroundColor Yellow
    Write-Host ""
}

# Step 4: Copy required DLLs from MSYS2 UCRT64 bin directory
Write-Host "Copying required DLLs from MSYS2 UCRT64..." -ForegroundColor Yellow
Write-Host ""

# Required DLL patterns (copy any matching files found)
$requiredDllPatterns = @(
    "libgcc_s_seh-1.dll",
    "libwinpthread-1.dll",
    "libgfortran-*.dll",
    "libquadmath-*.dll",
    "libgomp-1.dll",
    "libstdc++-6.dll",
    "libopenblas*.dll",
    "libfftw3-3.dll",
    "libfftw3_omp-3.dll"
)

$copiedDllCount = 0
$missingDlls = @()
$stagedDlls = @()

foreach ($pattern in $requiredDllPatterns) {
    $dllFiles = Get-ChildItem -Path $MingwBinDir -Filter $pattern -ErrorAction SilentlyContinue
    
    if ($dllFiles) {
        foreach ($dllFile in $dllFiles) {
            $destPath = Join-Path $DistDir $dllFile.Name
            try {
                Copy-Item -Path $dllFile.FullName -Destination $destPath -Force
                Write-Host "  Copied: $($dllFile.Name)" -ForegroundColor Green
                $stagedDlls += $dllFile.Name
                $copiedDllCount++
            } catch {
                Write-Host "  ERROR copying $($dllFile.Name): $_" -ForegroundColor Red
                $missingDlls += $dllFile.Name
            }
        }
    } else {
        Write-Host "  WARNING: Not found: $pattern" -ForegroundColor Yellow
        $missingDlls += $pattern
    }
}

Write-Host ""
Write-Host "Copied $copiedDllCount DLL(s)" -ForegroundColor Cyan
if ($missingDlls.Count -gt 0) {
    Write-Host "Missing patterns: $($missingDlls.Count)" -ForegroundColor Yellow
}
Write-Host ""

# Summary
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "Summary" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "  Executables copied: $copiedExeCount" -ForegroundColor Green
Write-Host "  DLLs copied: $copiedDllCount" -ForegroundColor Green
Write-Host ""
Write-Host "Destination (distribution) directory:" -ForegroundColor Cyan
Write-Host "  $DistDir" -ForegroundColor White
Write-Host ""

# List staged DLLs
if ($stagedDlls.Count -gt 0) {
    Write-Host "Staged DLLs:" -ForegroundColor Cyan
    foreach ($dll in $stagedDlls | Sort-Object) {
        Write-Host "  - $dll" -ForegroundColor Gray
    }
    Write-Host ""
}

# Verify distribution directory contents
Write-Host "Distribution directory contents:" -ForegroundColor Cyan
$distFiles = Get-ChildItem -Path $DistDir -File | Select-Object Name, Extension
$exeFiles = ($distFiles | Where-Object { $_.Extension -eq ".exe" }).Count
$dllFiles = ($distFiles | Where-Object { $_.Extension -eq ".dll" }).Count
Write-Host "  - $exeFiles executable(s)" -ForegroundColor Green
Write-Host "  - $dllFiles DLL(s)" -ForegroundColor Green

if ($exeFiles -eq 0) {
    Write-Host ""
    Write-Host "WARNING: No executables found in distribution directory!" -ForegroundColor Yellow
    Write-Host "Make sure QE has been built and installed to: $ArtifactsDir" -ForegroundColor Yellow
}

Write-Host ""

exit 0
