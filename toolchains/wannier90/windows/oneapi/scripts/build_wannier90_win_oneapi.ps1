# Wannier90 Windows Build Script for Intel oneAPI (CMake)
# This script builds Wannier90 using Intel oneAPI compilers (ifx/icx) via CMake
# Usage: .\scripts\build_wannier90_win_oneapi.ps1 [-NoMpi] [-W90SourceDir <path>]
#
# Prerequisites:
#   - Wannier90 source must be cloned from git (develop branch for CMake support)
#   - Run scripts/refresh_wannier90_source.ps1 -W90Version develop first

param(
    [switch]$NoMpi,                    # If specified, build serial Wannier90 without MPI
    [string]$W90SourceDir = "upstream/wannier90"  # Path to Wannier90 source directory (relative to repo root)
)

$ErrorActionPreference = 'Stop'

# Determine script location and repo root
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Resolve-Path (Join-Path $ScriptDir "..")

# Resolve Wannier90 source directory (can be relative to repo root or absolute)
if ([System.IO.Path]::IsPathRooted($W90SourceDir)) {
    $W90Source = Resolve-Path $W90SourceDir -ErrorAction Stop
} else {
    $W90Source = Resolve-Path (Join-Path $RepoRoot $W90SourceDir) -ErrorAction Stop
}

# Build directory is relative to Wannier90 source (not repo root)
$BuildDir = Join-Path $W90Source "build-win-oneapi"

Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "Wannier90 Windows oneAPI Build (CMake)" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Script directory: $ScriptDir"
Write-Host "Repository root: $RepoRoot"
Write-Host "Wannier90 source root: $W90Source"
Write-Host "Build directory: $BuildDir"
Write-Host ""

# Step 0: Verify Wannier90 source and CMakeLists.txt
Write-Host "Step 0: Verifying Wannier90 source..." -ForegroundColor Yellow

# Check if Wannier90 source directory exists
if (-not (Test-Path $W90Source)) {
    Write-Host ""
    Write-Host "ERROR: Wannier90 source directory not found: $W90Source" -ForegroundColor Red
    Write-Host ""
    Write-Host "Please run the refresh script first to clone Wannier90:" -ForegroundColor Yellow
    Write-Host "  .\scripts\refresh_wannier90_source.ps1 -W90Version develop" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "Or specify a different Wannier90 source directory with -W90SourceDir <path>" -ForegroundColor Yellow
    Write-Host ""
    exit 1
}

# Check for CMakeLists.txt (indicates this is a Wannier90 source directory with CMake support)
if (-not (Test-Path (Join-Path $W90Source "CMakeLists.txt"))) {
    Write-Host ""
    Write-Host "ERROR: Wannier90 source directory does not appear to have CMake support." -ForegroundColor Red
    Write-Host "  Missing: $W90Source\CMakeLists.txt" -ForegroundColor Red
    Write-Host ""
    Write-Host "CMake support is available in the develop branch." -ForegroundColor Yellow
    Write-Host "Please run: .\scripts\refresh_wannier90_source.ps1" -ForegroundColor Cyan
    Write-Host "  (The refresh script defaults to 'develop' branch which includes CMakeLists.txt)" -ForegroundColor Gray
    Write-Host ""
    exit 1
}

Write-Host "  Found Wannier90 source with CMake support" -ForegroundColor Green
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
# CMake is required for building Wannier90. We check PATH first, then common install locations.
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
# Wannier90 uses out-of-source builds, so we create a separate build directory
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
    "-S", "`"$($W90Source)`"",           # Source directory (Wannier90 root)
    "-B", "`"$($BuildDir)`"",            # Build directory
    "-G", "`"$cmakeGenerator`"",       # Generator (Ninja)
    "-DCMAKE_BUILD_TYPE=Release",        # Build type
    "-DCMAKE_Fortran_COMPILER=`"$ifxPath`"",    # Intel Fortran compiler (full path)
    "-DCMAKE_C_COMPILER=`"$icxPath`"",          # Intel C compiler (full path)
    "-DCMAKE_Fortran_COMPILER_FORCED=TRUE", # Skip Fortran compiler test (known CMake issue with ifx)
    "-DWANNIER90_MPI=$mpiFlag",         # MPI support (ON/OFF)
    "-DWANNIER90_SHARED_LIBS=OFF",      # Build static library
    "-DWANNIER90_INSTALL=OFF",           # Don't install (just build)
    "-DWANNIER90_TEST=OFF"              # Disable tests for initial build
) -join " "

# Build CMake build command
# Build all targets (wannier90.x, postw90.x, and library)
# Use --verbose to show compilation commands in real-time
$cmakeBuild = @(
    $cmakeExe,
    "--build", "`"$($BuildDir)`"",       # Build directory
    "--config", "Release",               # Configuration (for multi-config generators)
    "--parallel",                        # Use parallel build
    "--verbose"                          # Show compilation commands in real-time
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
Write-Host "Step 9: Locating built executables..." -ForegroundColor Yellow

# Wannier90 CMake build outputs executables with .x suffix
# Find executables in build directory
$exePatterns = @("wannier90.x", "postw90.x", "w90chk2chk.x", "w90spn2spn.x")
$foundExes = @()

foreach ($pattern in $exePatterns) {
    $exePath = Join-Path $BuildDir $pattern
    if (Test-Path $exePath) {
        $foundExes += $exePath
        Write-Host "  Found: $pattern" -ForegroundColor Green
    }
}

# Also check in bin subdirectory (some CMake configurations put executables there)
$binDir = Join-Path $BuildDir "bin"
if (Test-Path $binDir) {
    foreach ($pattern in $exePatterns) {
        $exePath = Join-Path $binDir $pattern
        if (Test-Path $exePath) {
            $foundExes += $exePath
            Write-Host "  Found: bin\$pattern" -ForegroundColor Green
        }
    }
}

if ($foundExes.Count -eq 0) {
    Write-Host "  WARNING: No executables found in expected locations" -ForegroundColor Yellow
    Write-Host "  Build directory: $BuildDir" -ForegroundColor Gray
} else {
    Write-Host "  Found $($foundExes.Count) executable(s)" -ForegroundColor Green
}
Write-Host ""

Write-Host "==========================================" -ForegroundColor Green
Write-Host "Build completed successfully!" -ForegroundColor Green
Write-Host "==========================================" -ForegroundColor Green
Write-Host ""
Write-Host "Executables should be in: $BuildDir" -ForegroundColor Cyan
if (Test-Path $binDir) {
    Write-Host "  or: $binDir" -ForegroundColor Cyan
}
Write-Host ""

exit 0

