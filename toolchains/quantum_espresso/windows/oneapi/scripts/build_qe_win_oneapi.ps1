# Quantum ESPRESSO Windows Build Script for Intel oneAPI
# This script builds QE using Intel oneAPI compilers (ifx/icx) via CMake
# Usage: .\scripts\build_qe_win_oneapi.ps1 [-NoMpi] [-QeSourceDir <path>]
#
# Prerequisites:
#   - QE source must be cloned from git (not a zip archive) with submodules initialized
#   - Run scripts/refresh_qe_source.ps1 first to set up upstream/qe/

param(
    [switch]$NoMpi,                    # If specified, build serial QE without MPI
    [string]$QeSourceDir = "upstream/qe"  # Path to QE source directory (relative to repo root)
)

$ErrorActionPreference = 'Stop'

# Determine script location and repo root
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Resolve-Path (Join-Path $ScriptDir "..")

# Resolve QE source directory (can be relative to repo root or absolute)
if ([System.IO.Path]::IsPathRooted($QeSourceDir)) {
    $QESource = Resolve-Path $QeSourceDir -ErrorAction Stop
} else {
    $QESource = Resolve-Path (Join-Path $RepoRoot $QeSourceDir) -ErrorAction Stop
}

# Build directory is relative to QE source (not repo root)
$BuildDir = Join-Path $QESource "build-win-oneapi"

Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "Quantum ESPRESSO Windows oneAPI Build" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Script directory: $ScriptDir"
Write-Host "Repository root: $RepoRoot"
Write-Host "QE source root: $QESource"
Write-Host "Build directory: $BuildDir"
Write-Host ""

# Step 0: Verify QE source and submodules
Write-Host "Step 0: Verifying QE source and submodules..." -ForegroundColor Yellow

# Check if QE source directory exists
if (-not (Test-Path $QESource)) {
    Write-Host ""
    Write-Host "ERROR: QE source directory not found: $QESource" -ForegroundColor Red
    Write-Host ""
    Write-Host "Please run the refresh script first to clone QE:" -ForegroundColor Yellow
    Write-Host "  .\scripts\refresh_qe_source.ps1" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "Or specify a different QE source directory with -QeSourceDir <path>" -ForegroundColor Yellow
    Write-Host ""
    exit 1
}

# Check for CMakeLists.txt (indicates this is a QE source directory)
if (-not (Test-Path (Join-Path $QESource "CMakeLists.txt"))) {
    Write-Host ""
    Write-Host "ERROR: QE source directory does not appear to be valid QE source." -ForegroundColor Red
    Write-Host "  Missing: $QESource\CMakeLists.txt" -ForegroundColor Red
    Write-Host ""
    Write-Host "This directory should be a git clone of QE with submodules initialized." -ForegroundColor Yellow
    Write-Host "Please run: .\scripts\refresh_qe_source.ps1" -ForegroundColor Cyan
    Write-Host ""
    exit 1
}

# Check for required submodules
$requiredSubmodules = @(
    @{Path = "external/wannier90"; CheckFile = "src/parameters.F90"},
    @{Path = "external/mbd"; CheckFile = "CMakeLists.txt"},
    @{Path = "external/devxlib"; CheckFile = "src/deviceXlib_mod.f90"}
)

$missingSubmodules = @()
foreach ($submodule in $requiredSubmodules) {
    $submodulePath = Join-Path $QESource $submodule.Path
    $checkFilePath = Join-Path $submodulePath $submodule.CheckFile
    
    if (-not (Test-Path $checkFilePath)) {
        $missingSubmodules += $submodule.Path
        Write-Host "  MISSING: $($submodule.Path) (expected file: $($submodule.CheckFile))" -ForegroundColor Red
    } else {
        Write-Host "  OK: $($submodule.Path)" -ForegroundColor Green
    }
}

if ($missingSubmodules.Count -gt 0) {
    Write-Host ""
    Write-Host "ERROR: Required submodules are missing or not properly initialized!" -ForegroundColor Red
    Write-Host ""
    Write-Host "Missing submodules: $($missingSubmodules -join ', ')" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "This usually happens when:" -ForegroundColor Yellow
    Write-Host "  1. QE was downloaded as a zip archive (not cloned from git)" -ForegroundColor Yellow
    Write-Host "  2. Submodules were not initialized during clone" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Solution:" -ForegroundColor Cyan
    Write-Host "  1. Run: .\scripts\refresh_qe_source.ps1" -ForegroundColor Cyan
    Write-Host "     This will clone QE from git with submodules properly initialized." -ForegroundColor Gray
    Write-Host ""
    Write-Host "  2. Or if you have a git clone, initialize submodules manually:" -ForegroundColor Cyan
    Write-Host "     cd $QESource" -ForegroundColor Gray
    Write-Host "     git submodule update --init --recursive" -ForegroundColor Gray
    Write-Host ""
    exit 1
}
Write-Host ""

# Step 1: Discover oneAPI installation
# We need to find setvars.bat to initialize the oneAPI environment
Write-Host "Step 1: Discovering Intel oneAPI..." -ForegroundColor Yellow

# Check if ONEAPI_ROOT is set in environment, otherwise use default location
if ($env:ONEAPI_ROOT) {
    $oneApiRoot = $env:ONEAPI_ROOT
    Write-Host "  Using ONEAPI_ROOT from environment: $oneApiRoot"
} else {
    $oneApiRoot = "C:\Program Files (x86)\Intel\oneAPI"
    Write-Host "  Trying default path: $oneApiRoot"
}

if (-not (Test-Path $oneApiRoot)) {
    throw "ERROR: Intel oneAPI not found at: $oneApiRoot. Please install Intel oneAPI HPC Toolkit or set ONEAPI_ROOT environment variable."
}

# Locate setvars.bat which initializes oneAPI environment variables
$setvarsBat = Join-Path $oneApiRoot "setvars.bat"
if (-not (Test-Path $setvarsBat)) {
    throw "ERROR: setvars.bat not found at: $setvarsBat"
}

Write-Host "  Found setvars.bat: $setvarsBat" -ForegroundColor Green
Write-Host ""

# Ensure MKLROOT and CMAKE_PREFIX_PATH are set for FindLAPACK
Write-Host "Step 1b: Preparing MKL environment for CMake..." -ForegroundColor Yellow
if ($env:MKLROOT -and (Test-Path $env:MKLROOT)) {
    $mklRoot = $env:MKLROOT
    Write-Host "  Using MKLROOT from environment: $mklRoot"
} else {
    $mklRoot = Join-Path $oneApiRoot "mkl\latest"
    Write-Host "  Deriving MKLROOT from ONEAPI_ROOT: $mklRoot"
}
if (-not (Test-Path $mklRoot)) {
    throw "ERROR: MKL not found at $mklRoot. Ensure oneAPI MKL is installed."
}
$env:MKLROOT = $mklRoot
$mklPrefix = Join-Path $mklRoot "lib\cmake\mkl"
$env:CMAKE_PREFIX_PATH = "$mklPrefix;$env:CMAKE_PREFIX_PATH"
Write-Host "  MKLROOT set to: $env:MKLROOT"
Write-Host "  CMAKE_PREFIX_PATH includes: $mklPrefix"
Write-Host ""

# Step 2: Determine build configuration
# Set MPI flag based on -NoMpi parameter
$mpiFlag = if ($NoMpi) { "OFF" } else { "ON" }
Write-Host "Step 2: Configuration" -ForegroundColor Yellow
Write-Host "  MPI: $mpiFlag"
Write-Host ""

# Step 3: Check for Visual Studio Build Tools (required for linker)
# Intel compilers on Windows need MSVC linker from Visual Studio Build Tools
Write-Host "Step 3: Checking for Visual Studio Build Tools..." -ForegroundColor Yellow
$vsFound = $false
$vsBasePath = $null
$vsPaths = @(
    "C:\Program Files\Microsoft Visual Studio\18",      # VS 2026
    "C:\Program Files\Microsoft Visual Studio\2022",     # VS 2022
    "C:\Program Files (x86)\Microsoft Visual Studio\2022",
    "C:\Program Files\Microsoft Visual Studio\2019",
    "C:\Program Files (x86)\Microsoft Visual Studio\2019"
)
foreach ($vsPath in $vsPaths) {
    if (Test-Path $vsPath) {
        $vsFound = $true
        $vsBasePath = $vsPath
        Write-Host "  Found Visual Studio at: $vsPath" -ForegroundColor Green
        break
    }
}
if (-not $vsFound) {
    Write-Host "  WARNING: Visual Studio Build Tools not found in standard locations" -ForegroundColor Yellow
    Write-Host "  The build may fail if the linker is not available." -ForegroundColor Yellow
    Write-Host "  Install from: https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2022" -ForegroundColor Cyan
}
Write-Host ""

# Step 4: Locate CMake executable
# CMake is required for building QE. We check PATH first, then common install locations.
Write-Host "Step 4: Checking for CMake..." -ForegroundColor Yellow
$cmakePath = $null

# Try to find cmake in PATH first
try {
    $cmakePath = (Get-Command cmake -ErrorAction Stop).Source
    Write-Host "  Found CMake in PATH: $cmakePath" -ForegroundColor Green
} catch {
    # If not in PATH, check common installation locations
    $commonPaths = @(
        "C:\Program Files\CMake\bin\cmake.exe",
        "C:\Program Files (x86)\CMake\bin\cmake.exe",
        "$env:ProgramFiles\CMake\bin\cmake.exe"
    )
    
    foreach ($p in $commonPaths) {
        if (Test-Path $p) {
            $cmakePath = $p
            # Add to PATH for this session so cmd.exe can find it
            $env:Path = "$(Split-Path $p);$env:Path"
            Write-Host "  Found CMake at: $cmakePath" -ForegroundColor Green
            Write-Host "  Added to PATH for this session" -ForegroundColor Yellow
            break
        }
    }
    
    if (-not $cmakePath) {
        throw "ERROR: CMake not found. Please install CMake and ensure it's in your PATH, or install it to a standard location (C:\Program Files\CMake)."
    }
}
Write-Host ""

# Step 5: Create build directory
# QE uses out-of-source builds, so we create a separate build directory
# If CMakeCache.txt exists from a previous generator, remove it to avoid conflicts
Write-Host "Step 5: Preparing build directory..." -ForegroundColor Yellow
if (Test-Path $BuildDir) {
    $cacheFile = Join-Path $BuildDir "CMakeCache.txt"
    if (Test-Path $cacheFile) {
        Write-Host "  Removing existing CMakeCache.txt (may be from different generator)" -ForegroundColor Yellow
        Remove-Item $cacheFile -Force -ErrorAction SilentlyContinue
        # Also remove CMakeFiles directory if it exists
        $cmakeFilesDir = Join-Path $BuildDir "CMakeFiles"
        if (Test-Path $cmakeFilesDir) {
            Remove-Item $cmakeFilesDir -Recurse -Force -ErrorAction SilentlyContinue
        }
    }
}
New-Item -ItemType Directory -Force -Path $BuildDir | Out-Null
Write-Host "  Created/verified: $BuildDir" -ForegroundColor Green
Write-Host ""

# Step 6: Determine CMake generator
# Ninja is the preferred generator for fast builds. We require it for this script.
Write-Host "Step 6: Determining CMake generator..." -ForegroundColor Yellow
$cmakeGenerator = $null

# Check if Ninja is available (required for this script)
try {
    $null = Get-Command ninja -ErrorAction Stop
    $cmakeGenerator = "Ninja"
    Write-Host "  Using Ninja generator" -ForegroundColor Green
} catch {
    Write-Host "  Ninja not found" -ForegroundColor Yellow
}

# If Ninja not found, provide clear installation instructions
if (-not $cmakeGenerator) {
    Write-Host ""
    Write-Host "ERROR: No suitable CMake generator found!" -ForegroundColor Red
    Write-Host ""
    Write-Host "Please install Ninja (required):" -ForegroundColor Yellow
    Write-Host "  1. Using Chocolatey: choco install ninja" -ForegroundColor Cyan
    Write-Host "  2. Or download from: https://github.com/ninja-build/ninja/releases" -ForegroundColor Cyan
    Write-Host "     Extract ninja.exe and add it to your PATH" -ForegroundColor Cyan
    Write-Host ""
    throw "No CMake generator available. Please install Ninja."
}
Write-Host ""

# Step 7: Build CMake commands
# Assemble the CMake configure and build commands as strings for cmd.exe
Write-Host "Step 7: Building CMake commands..." -ForegroundColor Yellow

# Use full path to cmake for cmd.exe (since PATH may not be set in cmd session)
$cmakeExe = "`"$cmakePath`""

# Find full paths to compilers (needed because PATH may not be set correctly in cmd.exe)
# Even though setvars.bat should set PATH, CMake may run before PATH is fully propagated
$ifxPath = $null
$icxPath = $null
$ifxSearch = Get-ChildItem "$oneApiRoot\compiler" -Recurse -Filter "ifx.exe" -ErrorAction SilentlyContinue | Select-Object -First 1
if ($ifxSearch) { 
    $ifxPath = $ifxSearch.FullName
    Write-Host "  Found ifx at: $ifxPath" -ForegroundColor Green
}
$icxSearch = Get-ChildItem "$oneApiRoot\compiler" -Recurse -Filter "icx.exe" -ErrorAction SilentlyContinue | Select-Object -First 1
if ($icxSearch) { 
    $icxPath = $icxSearch.FullName
    Write-Host "  Found icx at: $icxPath" -ForegroundColor Green
}

if (-not $ifxPath -or -not $icxPath) {
    throw "ERROR: Could not find ifx or icx compilers in oneAPI installation at: $oneApiRoot"
}
Write-Host ""

# Build CMake configure command
# Ninja is a single-config generator, so we use CMAKE_BUILD_TYPE
# Use full paths to compilers to ensure CMake can find them
$cmakeConfigure = @(
    $cmakeExe,
    "-S", "`"$($QESource)`"",           # Source directory (QE root)
    "-B", "`"$($BuildDir)`"",          # Build directory
    "-G", "`"$cmakeGenerator`"",       # Generator (Ninja)
    "-DCMAKE_BUILD_TYPE=Release",      # Build type
    "-DCMAKE_Fortran_COMPILER=`"$ifxPath`"",    # Intel Fortran compiler (full path)
    "-DCMAKE_C_COMPILER=`"$icxPath`"",          # Intel C compiler (full path)
    "-DCMAKE_Fortran_COMPILER_FORCED=TRUE", # Skip Fortran compiler test (known CMake issue with ifx)
    "-DQE_CPP=`"$icxPath`"",                    # Use icx as C preprocessor (QE requirement)
    "-DQE_ENABLE_MPI=$mpiFlag",        # MPI support (ON/OFF)
    "-DQE_LAPACK_INTERNAL=OFF",        # Use external MKL, not internal LAPACK
    "-DQE_ENABLE_OPENMP=OFF",          # Disable OpenMP for initial build
    "-DQE_ENABLE_SCALAPACK=OFF",       # Disable SCALAPACK for initial build
    "-DQE_ENABLE_HDF5=OFF",            # Disable HDF5 for initial build
    "-DQE_ENABLE_LIBXC=OFF",           # Disable LibXC for initial build
    "-DQE_ENABLE_ELPA=OFF",            # Disable ELPA for initial build
    "-DQE_ENABLE_TEST=OFF",            # Disable tests for initial build
    "-DQE_WANNIER90_INTERNAL=ON",      # Use internal Wannier90 (submodule exists)
    "-DQE_MBD_INTERNAL=ON",            # Use internal MBD (submodule exists)
    "-DQE_DEVICEXLIB_INTERNAL=ON",     # Use internal DeviceXlib (submodule exists)
    "-DQE_ENABLE_ENVIRON=NO",          # Disable Environ (avoids git submodule issues)
    "-DQE_FFTW_VENDOR=Internal"        # Use internal FFTW (MKL FFTW detection may be complex)
) -join " "

# Build CMake build command
# Build only the 'pw' target initially (can be extended later)
# Use --verbose to show compilation commands in real-time
$cmakeBuild = @(
    $cmakeExe,
    "--build", "`"$($BuildDir)`"",     # Build directory
    "--target", "pw",                  # Target to build (pw.x executable)
    "--config", "Release",             # Configuration (for multi-config generators)
    "--parallel",                       # Use parallel build
    "--verbose"                         # Show compilation commands in real-time
) -join " "

# Step 8: Execute build in single cmd.exe process
# IMPORTANT: All commands must run in the same cmd.exe session so that:
# 1. setvars.bat initializes oneAPI environment (ifx, icx, MKL in PATH)
# 2. CMake configure uses those compilers
# 3. CMake build uses those compilers
# If we ran cmake directly in PowerShell, it wouldn't have the oneAPI environment.
Write-Host "Step 8: Running build..." -ForegroundColor Yellow
Write-Host ""

# Find Visual Studio initialization script (needed for linker)
# VS 2022+ uses VsDevCmd.bat, older versions use vcvarsall.bat
$vsInitScript = $null
if ($vsFound -and $vsBasePath) {
    # Check all editions in the found VS installation
    $vsEditions = Get-ChildItem $vsBasePath -Directory -ErrorAction SilentlyContinue | Select-Object -ExpandProperty Name
    foreach ($edition in $vsEditions) {
        # Try VsDevCmd.bat first (VS 2022+, recommended)
        $vsDevCmdPath = Join-Path $vsBasePath "$edition\Common7\Tools\VsDevCmd.bat"
        if (Test-Path $vsDevCmdPath) {
            $vsInitScript = $vsDevCmdPath
            Write-Host "  Found VsDevCmd.bat: $vsInitScript" -ForegroundColor Green
            break
        }
        # Fall back to vcvarsall.bat (older VS versions or if VsDevCmd not found)
        $vcvarsPath = Join-Path $vsBasePath "$edition\VC\Auxiliary\Build\vcvarsall.bat"
        if (Test-Path $vcvarsPath) {
            $vsInitScript = $vcvarsPath
            Write-Host "  Found vcvarsall.bat: $vsInitScript" -ForegroundColor Green
            break
        }
    }
}

# Note: Commands are assembled in the try block below to show output separately
if ($vsInitScript) {
    Write-Host "  Using Visual Studio linker initialization" -ForegroundColor Green
} else {
    Write-Host "  Using setvars.bat VS initialization (may fail if VS not detected)" -ForegroundColor Yellow
}
Write-Host ""

# Execute in cmd.exe and capture exit code
# Temporarily disable error action to allow CMake's informational messages
$oldErrorAction = $ErrorActionPreference
$ErrorActionPreference = 'Continue'

try {
    # Run the command and stream output in real-time
    Write-Host "Starting CMake configuration and build..." -ForegroundColor Cyan
    Write-Host ""
    
    # Build initialization prefix (VS + oneAPI setvars)
    if ($vsInitScript) {
        if ($vsInitScript -like "*VsDevCmd.bat") {
            $vsInitPrefix = "`"$vsInitScript`" -arch=x64 -host_arch=x64 >nul 2>&1"
        } else {
            $vsInitPrefix = "`"$vsInitScript`" x64 >nul 2>&1"
        }
    } else {
        $vsInitPrefix = "`"$setvarsBat`" intel64 vs2022 >nul 2>&1"
    }
    $oneApiInit = "`"$setvarsBat`" intel64 >nul 2>&1"
    
    # Split configure and build to show output separately with verbose flags
    $configureCmd = "$vsInitPrefix && $oneApiInit && $cmakeConfigure"
    $buildCmd = "$vsInitPrefix && $oneApiInit && $cmakeBuild"
    
    Write-Host "=== CMake Configuration ===" -ForegroundColor Yellow
    & cmd.exe /c $configureCmd
    $configureExit = $LASTEXITCODE
    if ($configureExit -ne 0) {
        throw "CMake configuration failed with exit code $configureExit"
    }
    
    Write-Host ""
    Write-Host "=== CMake Build (verbose) ===" -ForegroundColor Yellow
    & cmd.exe /c $buildCmd
    $exitCode = $LASTEXITCODE
    
    # Check if build succeeded
    if ($exitCode -ne 0) {
        Write-Host ""
        Write-Host "ERROR: Build failed with exit code $exitCode" -ForegroundColor Red
        exit $exitCode
    }
} finally {
    # Restore error action preference
    $ErrorActionPreference = $oldErrorAction
}

Write-Host ""
Write-Host "==========================================" -ForegroundColor Green
Write-Host "Build completed successfully!" -ForegroundColor Green
Write-Host "==========================================" -ForegroundColor Green
Write-Host ""
Write-Host "Executables should be in: $BuildDir\bin" -ForegroundColor Cyan
Write-Host ""

exit 0

