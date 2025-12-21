# Quantum ESPRESSO Windows Build Script for Intel oneAPI
# This script builds QE using Intel oneAPI compilers (ifx/icx) via CMake
# Usage: .\scripts\build_qe_win_oneapi.ps1 [-NoMpi] [-QeSourceDir <path>] [-FftwVendor <Internal|Intel_FFTW3|Intel_DFTI|AUTO>]
#
# Prerequisites:
#   - QE source must be cloned from git (not a zip archive) with submodules initialized
#   - Run scripts/refresh_qe_source.ps1 first to set up upstream/qe/
#
# FFT Backend Options:
#   -FftwVendor Intel_FFTW3 (default): Use external Intel MKL FFTW3 wrapper
#   -FftwVendor Internal: Use QE's internal FFT implementation
#   -FftwVendor Intel_DFTI: Use Intel MKL DFTI interface
#   -FftwVendor AUTO: Let CMake auto-detect (not recommended)

param(
    [switch]$NoMpi,                    # If specified, build serial QE without MPI
    [string]$QeSourceDir = "upstream/qe",  # Path to QE source directory (relative to repo root)
    [ValidateSet('Internal', 'Intel_FFTW3', 'Intel_DFTI', 'AUTO')]
    [string]$FftwVendor = "Intel_DFTI",  # FFT backend vendor (default: Intel_DFTI)
    [ValidateSet('pw', 'all')]
    [string]$MakeTarget = "pw"  # Build target: 'pw' (pw only, default) or 'all' (all executables)
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
# Use separate build directory for MPI builds to avoid cache conflicts
$buildDirSuffix = if ($NoMpi) { "build-win-oneapi" } else { "build-win-oneapi-msmpi" }
$BuildDir = Join-Path $QESource $buildDirSuffix

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
Write-Host "  FFT Vendor: $FftwVendor"
Write-Host "  Build Target: $MakeTarget" -ForegroundColor Cyan
if ($MakeTarget -eq "all") {
    Write-Host "  Mode: Building ALL QE executables" -ForegroundColor Green
} else {
    Write-Host "  Mode: Building pw only (default)" -ForegroundColor Green
}
Write-Host ""

# Step 2b: Validate MS-MPI installation when MPI is enabled
$msMpiPaths = @{}
if ($mpiFlag -eq "ON") {
    Write-Host "Step 2b: Validating MS-MPI installation..." -ForegroundColor Yellow
    
    # Define MS-MPI SDK paths (no guessing beyond defaults)
    # Assume 64-bit build (VS Hostx64/x64 + ifx x64), so always choose Include\x64
    $msmpiSdk = "C:\Program Files (x86)\Microsoft SDKs\MPI"
    $msmpiIncBase = Join-Path $msmpiSdk "Include"
    $msmpiIncArch = Join-Path $msmpiSdk "Include\x64"
    $msmpiLibDir = Join-Path $msmpiSdk "Lib\x64"
    
    # Store paths for later use
    $msMpiPaths['MPIEXEC'] = "C:\Program Files\Microsoft MPI\Bin\mpiexec.exe"
    $msMpiPaths['SDK_ROOT'] = $msmpiSdk
    $msMpiPaths['SDK_INCLUDE_BASE'] = $msmpiIncBase
    $msMpiPaths['SDK_INCLUDE_ARCH'] = $msmpiIncArch
    $msMpiPaths['SDK_LIB_DIR'] = $msmpiLibDir
    $msMpiPaths['MPIF_H'] = Join-Path $msmpiIncBase "mpif.h"
    $msMpiPaths['MPIFPTR_H'] = Join-Path $msmpiIncArch "mpifptr.h"
    
    # Define MS-MPI library (only msmpi.lib needed for validation)
    # The FindMPI shim will handle all library linking
    $msmpiLib = Join-Path $msmpiLibDir "msmpi.lib"
    $msMpiPaths['SDK_LIB'] = $msmpiLib
    
    # Validate required files exist
    $missingPaths = @()
    
    # Check mpif.h in base Include directory
    if (-not (Test-Path $msMpiPaths['MPIF_H'])) {
        Write-Host "  MISSING: mpif.h = $($msMpiPaths['MPIF_H'])" -ForegroundColor Red
        $missingPaths += 'MPIF_H'
    } else {
        Write-Host "  Found: mpif.h = $($msMpiPaths['MPIF_H'])" -ForegroundColor Green
    }
    
    # Check mpifptr.h in architecture-specific Include\x64 directory
    if (-not (Test-Path $msMpiPaths['MPIFPTR_H'])) {
        Write-Host "  MISSING: mpifptr.h = $($msMpiPaths['MPIFPTR_H'])" -ForegroundColor Red
        $missingPaths += 'MPIFPTR_H'
    } else {
        Write-Host "  Found: mpifptr.h = $($msMpiPaths['MPIFPTR_H'])" -ForegroundColor Green
    }
    
    # Validate msmpi.lib exists (FindMPI shim will validate msmpifec.lib)
    if (-not (Test-Path $msmpiLib)) {
        Write-Host "  MISSING: msmpi.lib = $msmpiLib" -ForegroundColor Red
        $missingPaths += 'msmpi.lib'
    } else {
        Write-Host "  Found: msmpi.lib = $msmpiLib" -ForegroundColor Green
    }
    
    # Also validate directory paths
    foreach ($key in @('MPIEXEC', 'SDK_ROOT', 'SDK_INCLUDE_BASE', 'SDK_INCLUDE_ARCH', 'SDK_LIB_DIR')) {
        $path = $msMpiPaths[$key]
        if (Test-Path $path) {
            Write-Host "  Found: $key = $path" -ForegroundColor Green
        } else {
            Write-Host "  MISSING: $key = $path" -ForegroundColor Red
            $missingPaths += $key
        }
    }
    
    if ($missingPaths.Count -gt 0) {
        Write-Host ""
        Write-Host "ERROR: MS-MPI SDK incomplete; need mpif.h + mpifptr.h + msmpi.lib" -ForegroundColor Red
        Write-Host ""
        Write-Host "Missing components:" -ForegroundColor Yellow
        foreach ($key in $missingPaths) {
            if ($key -match '\.lib$') {
                Write-Host "  - $key" -ForegroundColor Yellow
            } else {
                Write-Host "  - $key : $($msMpiPaths[$key])" -ForegroundColor Yellow
            }
        }
        Write-Host ""
        Write-Host "Download from: https://www.microsoft.com/en-us/download/details.aspx?id=57467" -ForegroundColor Cyan
        Write-Host ""
        exit 1
    }
    
    # Workaround C: Copy mpifptr.h into QE include directory for MS-MPI Fortran builds
    Write-Host "Workaround C: copying mpifptr.h into QE include directory..." -ForegroundColor Yellow
    $src = Join-Path $msmpiSdk "Include\x64\mpifptr.h"
    $qeIncludeDir = Join-Path $QESource "include"
    $dst = Join-Path $qeIncludeDir "mpifptr.h"
    
    # Verify source exists
    if (-not (Test-Path $src)) {
        throw "ERROR: Workaround C failed - source file not found: $src"
    }
    
    # Ensure QE include directory exists
    if (-not (Test-Path $qeIncludeDir)) {
        throw "ERROR: Workaround C failed - QE include directory does not exist: $qeIncludeDir (check QE source root: $QESource)"
    }
    
    # Copy the file
    Copy-Item -Path $src -Destination $dst -Force -ErrorAction Stop
    
    # Verify destination exists after copy
    if (-not (Test-Path $dst)) {
        throw "ERROR: Workaround C failed - destination file was not created: $dst"
    }
    
    # Get file info and print success message
    $info = Get-Item $dst
    Write-Host "Workaround C: copied mpifptr.h into QE include directory" -ForegroundColor Green
    Write-Host "  Source: $src" -ForegroundColor White
    Write-Host "  Destination: $dst" -ForegroundColor White
    Write-Host "  File size: $($info.Length) bytes, LastWriteTime: $($info.LastWriteTime)" -ForegroundColor White
    Write-Host ""
}

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

# Step 6b: Resolve MKL library paths (needed for BLAS/LAPACK and VendorFFTW)
Write-Host "Step 6b: Resolving MKL library paths..." -ForegroundColor Yellow

# Resolve MKL libraries (always needed since we use external MKL for LAPACK)
$mklRtLib = Join-Path $mklRoot "lib\mkl_rt.lib"
if (-not (Test-Path $mklRtLib)) {
    Write-Host ""
    Write-Host "ERROR: MKL library not found!" -ForegroundColor Red
    Write-Host "  Missing: $mklRtLib" -ForegroundColor Red
    Write-Host ""
    exit 1
}

# Resolve libiomp5md.lib
$libiomp5mdPath = $null
$latestCompilerPath = Join-Path $oneApiRoot "compiler\latest\lib\libiomp5md.lib"
if (Test-Path $latestCompilerPath) {
    $libiomp5mdPath = $latestCompilerPath
} else {
    # Search for all compiler versions and pick the newest
    $compilerVersions = Get-ChildItem "$oneApiRoot\compiler" -Directory -ErrorAction SilentlyContinue | 
        Where-Object { $_.Name -match '^\d+\.\d+$' } |
        Sort-Object { [Version]$_.Name } -Descending
    
    foreach ($version in $compilerVersions) {
        $candidate = Join-Path $version.FullName "lib\libiomp5md.lib"
        if (Test-Path $candidate) {
            $libiomp5mdPath = $candidate
            break
        }
    }
}

if (-not $libiomp5mdPath) {
    Write-Host ""
    Write-Host "ERROR: libiomp5md.lib not found!" -ForegroundColor Red
    Write-Host "  Searched in: $oneApiRoot\compiler\*\lib\libiomp5md.lib" -ForegroundColor Red
    Write-Host ""
    Write-Host "MKL threading typically requires libiomp5md.lib." -ForegroundColor Yellow
    Write-Host ""
    exit 1
}

# Set BLAS/LAPACK libraries to full absolute paths (needed for FindVendorFFTW)
$blasLapackLibraries = "$mklRtLib;$libiomp5mdPath"
Write-Host "  BLAS/LAPACK libraries: $blasLapackLibraries" -ForegroundColor Green
Write-Host "  mkl_rt.lib: $mklRtLib" -ForegroundColor Green
Write-Host "  libiomp5md.lib: $libiomp5mdPath" -ForegroundColor Green
Write-Host ""

# Step 6c: Validate FFT vendor and resolve paths
Write-Host "Step 6c: Validating FFT vendor configuration..." -ForegroundColor Yellow
Write-Host "  FFT vendor: $FftwVendor" -ForegroundColor Cyan

$fftVendorIncludeDirs = $null
$fftVendorLibraries = $null

if ($FftwVendor -eq "Intel_FFTW3" -or $FftwVendor -eq "Intel_DFTI") {
    # Validate required paths for both Intel_FFTW3 and Intel_DFTI
    $fftw3Header = Join-Path $mklRoot "include\fftw\fftw3.h"
    
    if (-not (Test-Path $fftw3Header)) {
        Write-Host ""
        Write-Host "ERROR: $FftwVendor selected but required header not found!" -ForegroundColor Red
        Write-Host "  Missing: $fftw3Header" -ForegroundColor Red
        Write-Host ""
        Write-Host "Either install MKL FFTW headers/libs or rerun with -FftwVendor Internal" -ForegroundColor Yellow
        Write-Host ""
        exit 1
    }
    
    # Set paths for CMake (same for both Intel_FFTW3 and Intel_DFTI)
    $fftVendorIncludeDirs = Join-Path $mklRoot "include\fftw"
    $fftVendorLibraries = "$mklRtLib;$libiomp5mdPath"
    
    Write-Host "  VendorFFTW_INCLUDE_DIRS: $fftVendorIncludeDirs" -ForegroundColor Green
    Write-Host "  VendorFFTW_LIBRARIES: $fftVendorLibraries" -ForegroundColor Green
} elseif ($FftwVendor -eq "Internal") {
    Write-Host "  Using QE internal FFT implementation" -ForegroundColor Green
} elseif ($FftwVendor -eq "AUTO") {
    Write-Host "  Letting CMake auto-detect FFT backend" -ForegroundColor Yellow
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
    "-DCMAKE_Fortran_FLAGS_RELEASE=`"-O2 -DNDEBUG -heap-arrays`"",  # Release: -O2 optimization, heap arrays for stack overflow fix
    "-DCMAKE_C_FLAGS_RELEASE=`"-O2 -DNDEBUG`"",  # Release: -O2 optimization for C code
    "-DCMAKE_Fortran_FLAGS_DEBUG=`"-heap-arrays -traceback`"",  # Use heap arrays + traceback for Debug
    "-DCMAKE_EXE_LINKER_FLAGS=`"-Qoption,link,/STACK:134217728`"",  # Increase stack size to 128MB (fix stack overflow) - use Intel pass-through to avoid ifx warning
    "-DCMAKE_SHARED_LINKER_FLAGS=`"-Qoption,link,/STACK:134217728`"",  # Also apply to shared modules
    "-DCMAKE_Fortran_COMPILER_FORCED=TRUE", # Skip Fortran compiler test (known CMake issue with ifx)
    "-DQE_CPP=`"$icxPath`"",                    # Use icx as C preprocessor (QE requirement)
    "-DQE_ENABLE_MPI=$mpiFlag",        # MPI support (ON/OFF)
    "-DQE_LAPACK_INTERNAL=OFF",        # Use external MKL, not internal LAPACK
    "-DQE_ENABLE_OPENMP=ON",           # Enable OpenMP
    "-DQE_ENABLE_SCALAPACK=OFF",       # Disable SCALAPACK for initial build
    "-DQE_ENABLE_HDF5=OFF",            # Disable HDF5 for initial build
    "-DQE_ENABLE_LIBXC=ON",            # Enable LibXC for production
    "-DQE_ENABLE_ELPA=OFF",            # Disable ELPA for initial build
    "-DQE_ENABLE_TEST=OFF",            # Disable tests for initial build
    "-DQE_WANNIER90_INTERNAL=ON",      # Use internal Wannier90 (submodule exists)
    "-DQE_MBD_INTERNAL=ON",            # Use internal MBD (submodule exists)
    "-DQE_DEVICEXLIB_INTERNAL=ON",     # Use internal DeviceXlib (submodule exists)
    "-DQE_ENABLE_ENVIRON=NO",          # Disable Environ (avoids git submodule issues)
    "-DQE_FFTW_VENDOR=$FftwVendor"     # FFT vendor (user-selected)
)

# Force BLAS_LIBRARIES and LAPACK_LIBRARIES to full absolute paths
# This ensures FindVendorFFTW can build VendorFFTW_LIBRARIES from them correctly
$cmakeConfigure += "-DBLAS_LIBRARIES=`"$blasLapackLibraries`""
$cmakeConfigure += "-DLAPACK_LIBRARIES=`"$blasLapackLibraries`""

# Add VendorFFTW hints for both Intel_FFTW3 and Intel_DFTI
# Note: FindVendorFFTW module builds VendorFFTW_LIBRARIES from BLAS_LIBRARIES and LAPACK_LIBRARIES
# when it finds MKL in the search paths. However, we set VendorFFTW_LIBRARIES as required
# by the user specification. The module may override this, but setting it ensures the
# variables are available if the module needs them.
if ($FftwVendor -eq "Intel_FFTW3" -or $FftwVendor -eq "Intel_DFTI") {
    $cmakeConfigure += "-DVendorFFTW_ID=Intel"
    $cmakeConfigure += "-DVendorFFTW_INCLUDE_DIRS=`"$fftVendorIncludeDirs`""
    $cmakeConfigure += "-DVendorFFTW_LIBRARIES=`"$fftVendorLibraries`""
}

# Add MS-MPI configuration when MPI is enabled
if ($mpiFlag -eq "ON") {
    # Force QE to use mpif.h (MS-MPI typically does not provide Fortran mpi module)
    $cmakeConfigure += "-DQE_ENABLE_MPI_MODULE=OFF"
    
    # Add CMAKE_MODULE_PATH to use our custom FindMPI.cmake shim
    # This bypasses CMake's built-in FindMPI to avoid MPI_Fortran_WORKS try_compile tests
    $cmakeModulePath = Join-Path $ScriptDir "cmake"
    $cmakeModulePathCmake = $cmakeModulePath -replace '\\', '/'
    $cmakeConfigure += "-DCMAKE_MODULE_PATH=`"$cmakeModulePathCmake`""
    
    # Convert Windows paths to CMake-friendly format (forward slashes)
    $msMpiExec = $msMpiPaths['MPIEXEC'] -replace '\\', '/'
    $msMpiHome = $msMpiPaths['SDK_ROOT'] -replace '\\', '/'
    $msMpiIncBase = $msMpiPaths['SDK_INCLUDE_BASE'] -replace '\\', '/'
    $msMpiIncArch = $msMpiPaths['SDK_INCLUDE_ARCH'] -replace '\\', '/'
    
    # Critical: prevent FindMPI from selecting Intel mpiifx.bat by forcing the MPI Fortran compiler to be the real ifx
    $ifxPathCmake = $ifxPath -replace '\\', '/'
    
    # Pass both include directories (base and x64) - semicolon-separated for CMake
    $msMpiIncDirs = "$msMpiIncBase;$msMpiIncArch"
    
    # Force FindMPI to use MSMPI and skip compiler wrappers
    # Note: The FindMPI shim handles all library linking, so we don't pass MPI_*_LIBRARIES here
    $cmakeConfigure += "-DMPI_GUESS_LIBRARY_NAME=MSMPI"
    $cmakeConfigure += "-DMPI_SKIP_COMPILER_WRAPPER=TRUE"
    $cmakeConfigure += "-DMPIEXEC_EXECUTABLE=`"$msMpiExec`""
    $cmakeConfigure += "-DMPI_HOME=`"$msMpiHome`""
    $cmakeConfigure += "-DMPI_Fortran_COMPILER=`"$ifxPathCmake`""
    $cmakeConfigure += "-DMPI_Fortran_INCLUDE_DIRS=`"$msMpiIncDirs`""
}

# Convert to single string for cmd.exe
$cmakeConfigure = $cmakeConfigure -join " "

# Build CMake build command
# Build target based on MakeTarget parameter (pw or all)
# Use --verbose to show compilation commands in real-time
$buildTarget = if ($MakeTarget -eq "all") { "all" } else { "pw" }
Write-Host "Build target: $buildTarget" -ForegroundColor Cyan

$cmakeBuild = @(
    $cmakeExe,
    "--build", "`"$($BuildDir)`"",     # Build directory
    "--target", $buildTarget,          # Target to build (pw.x executable or all)
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
    
    # Set MS-MPI environment variables before CMake configure (only when MPI is enabled)
    if ($mpiFlag -eq "ON") {
        $env:MSMPI_INC = $msMpiPaths['SDK_INCLUDE_BASE']
        $env:MSMPI_LIB64 = $msMpiPaths['SDK_LIB_DIR']
        Write-Host "=== MS-MPI Environment Variables ===" -ForegroundColor Cyan
        Write-Host "MSMPI_INC: $env:MSMPI_INC" -ForegroundColor White
        Write-Host "MSMPI_LIB64: $env:MSMPI_LIB64" -ForegroundColor White
        Write-Host ""
        Write-Host "MS-MPI Fortran include dirs: $($msMpiPaths['SDK_INCLUDE_BASE']) ; $($msMpiPaths['SDK_INCLUDE_ARCH'])" -ForegroundColor Cyan
        Write-Host ""
    }
    
    # Log library configuration before running CMake
    Write-Host "=== Library Configuration ===" -ForegroundColor Cyan
    Write-Host "BLAS_LIBRARIES: $blasLapackLibraries" -ForegroundColor White
    Write-Host "LAPACK_LIBRARIES: $blasLapackLibraries" -ForegroundColor White
    Write-Host "FFTW vendor: $FftwVendor" -ForegroundColor White
    if ($fftVendorIncludeDirs) {
        Write-Host "VendorFFTW_INCLUDE_DIRS: $fftVendorIncludeDirs" -ForegroundColor White
    }
    if ($fftVendorLibraries) {
        Write-Host "VendorFFTW_LIBRARIES: $fftVendorLibraries" -ForegroundColor White
    }
    Write-Host ""
    
    Write-Host "=== CMake Configuration ===" -ForegroundColor Yellow
    # Capture configure output to verify MS-MPI shim was used
    # Use try-catch to suppress PowerShell's false error interpretation of cmd.exe output
    try {
        $configureOutput = & cmd.exe /c $configureCmd 2>&1 | Out-String
    } catch {
        # PowerShell may interpret some cmd.exe output as errors, but configure still succeeds
        # Re-run without error handling to get actual output
        $configureOutput = cmd.exe /c $configureCmd 2>&1 | Out-String
    }
    Write-Host $configureOutput
    $configureExit = $LASTEXITCODE
    if ($configureExit -ne 0) {
        Write-Host ""
        Write-Host "ERROR: CMake configuration failed with exit code $configureExit" -ForegroundColor Red
        if ($FftwVendor -eq "Intel_FFTW3" -or $FftwVendor -eq "Intel_DFTI") {
            Write-Host ""
            Write-Host "$FftwVendor selected but VendorFFTW not found." -ForegroundColor Yellow
            Write-Host "Either install MKL FFTW headers/libs or rerun with -FftwVendor Internal" -ForegroundColor Yellow
            Write-Host ""
        }
        exit $configureExit
    }
    
    # Verify MS-MPI shim was used (only when MPI is enabled)
    if ($mpiFlag -eq "ON") {
        Write-Host ""
        Write-Host "=== Verifying MS-MPI shim usage ===" -ForegroundColor Cyan
        if ($configureOutput -match "Using MS-MPI shim") {
            Write-Host "  OK: Found 'Using MS-MPI shim' in CMake output" -ForegroundColor Green
        } else {
            Write-Host "  WARNING: 'Using MS-MPI shim' not found in CMake output" -ForegroundColor Yellow
            Write-Host "  The custom FindMPI.cmake shim may not have been used" -ForegroundColor Yellow
        }
        
        # Check CMakeCache.txt does NOT mention impi / mpiifx
        $cmakeCacheFile = Join-Path $BuildDir "CMakeCache.txt"
        if (Test-Path $cmakeCacheFile) {
            $impiMatches = Select-String -Path $cmakeCacheFile -Pattern 'impi' -CaseSensitive:$false
            $mpiifxMatches = Select-String -Path $cmakeCacheFile -Pattern 'mpiifx' -CaseSensitive:$false
            
            if ($impiMatches) {
                Write-Host "  WARNING: Found 'impi' references in CMakeCache.txt (Intel MPI detected)" -ForegroundColor Yellow
                $impiMatches | ForEach-Object { Write-Host "    $($_.Line)" -ForegroundColor Gray }
            } else {
                Write-Host "  OK: No 'impi' references in CMakeCache.txt" -ForegroundColor Green
            }
            
            if ($mpiifxMatches) {
                Write-Host "  WARNING: Found 'mpiifx' references in CMakeCache.txt (Intel MPI wrapper detected)" -ForegroundColor Yellow
                $mpiifxMatches | ForEach-Object { Write-Host "    $($_.Line)" -ForegroundColor Gray }
            } else {
                Write-Host "  OK: No 'mpiifx' references in CMakeCache.txt" -ForegroundColor Green
            }
        } else {
            Write-Host "  WARNING: CMakeCache.txt not found for verification" -ForegroundColor Yellow
        }
        Write-Host ""
    }
    
    # Validate resolved library paths from CMakeCache.txt
    Write-Host ""
    Write-Host "=== Validating CMake Cache Variables ===" -ForegroundColor Cyan
    $cmakeCacheFile = Join-Path $BuildDir "CMakeCache.txt"
    if (Test-Path $cmakeCacheFile) {
        # Extract BLAS_LIBRARIES
        $blasLine = Select-String -Path $cmakeCacheFile -Pattern '^BLAS_LIBRARIES:STRING=' | Select-Object -First 1
        if ($blasLine) {
            $blasValue = $blasLine.Line -replace '^BLAS_LIBRARIES:STRING=', ''
            Write-Host "BLAS_LIBRARIES: $blasValue" -ForegroundColor White
        } else {
            Write-Host "BLAS_LIBRARIES: (not found in cache)" -ForegroundColor Yellow
        }
        
        # Extract LAPACK_LIBRARIES
        $lapackLine = Select-String -Path $cmakeCacheFile -Pattern '^LAPACK_LIBRARIES:STRING=' | Select-Object -First 1
        if ($lapackLine) {
            $lapackValue = $lapackLine.Line -replace '^LAPACK_LIBRARIES:STRING=', ''
            Write-Host "LAPACK_LIBRARIES: $lapackValue" -ForegroundColor White
        } else {
            Write-Host "LAPACK_LIBRARIES: (not found in cache)" -ForegroundColor Yellow
        }
        
        # Extract VendorFFTW_LIBRARIES
        $vendorLine = Select-String -Path $cmakeCacheFile -Pattern '^VendorFFTW_LIBRARIES:STRING=' | Select-Object -First 1
        if ($vendorLine) {
            $vendorValue = $vendorLine.Line -replace '^VendorFFTW_LIBRARIES:STRING=', ''
            Write-Host "VendorFFTW_LIBRARIES: $vendorValue" -ForegroundColor White
        } else {
            Write-Host "VendorFFTW_LIBRARIES: (not found in cache)" -ForegroundColor Yellow
        }
        
        # Extract CMAKE_Fortran_FLAGS_RELEASE
        $fortranFlagsLine = Select-String -Path $cmakeCacheFile -Pattern '^CMAKE_Fortran_FLAGS_RELEASE:STRING=' | Select-Object -First 1
        if ($fortranFlagsLine) {
            $fortranFlagsValue = $fortranFlagsLine.Line -replace '^CMAKE_Fortran_FLAGS_RELEASE:STRING=', ''
            Write-Host "CMAKE_Fortran_FLAGS_RELEASE: $fortranFlagsValue" -ForegroundColor White
        } else {
            Write-Host "CMAKE_Fortran_FLAGS_RELEASE: (not found in cache)" -ForegroundColor Yellow
        }
        
        # Extract CMAKE_C_FLAGS_RELEASE
        $cFlagsLine = Select-String -Path $cmakeCacheFile -Pattern '^CMAKE_C_FLAGS_RELEASE:STRING=' | Select-Object -First 1
        if ($cFlagsLine) {
            $cFlagsValue = $cFlagsLine.Line -replace '^CMAKE_C_FLAGS_RELEASE:STRING=', ''
            Write-Host "CMAKE_C_FLAGS_RELEASE: $cFlagsValue" -ForegroundColor White
        } else {
            Write-Host "CMAKE_C_FLAGS_RELEASE: (not found in cache)" -ForegroundColor Yellow
        }
        
        # Extract MPI-related cache values (post-configure sanity check)
        if ($mpiFlag -eq "ON") {
            # Extract MPI_GUESS_LIBRARY_NAME
            $mpiGuessLibLine = Select-String -Path $cmakeCacheFile -Pattern '^MPI_GUESS_LIBRARY_NAME:STRING=' | Select-Object -First 1
            if ($mpiGuessLibLine) {
                $mpiGuessLibValue = $mpiGuessLibLine.Line -replace '^MPI_GUESS_LIBRARY_NAME:STRING=', ''
                Write-Host "MPI_GUESS_LIBRARY_NAME: $mpiGuessLibValue" -ForegroundColor White
            } else {
                Write-Host "MPI_GUESS_LIBRARY_NAME: (not found in cache)" -ForegroundColor Yellow
            }
            
            # Extract MPI_SKIP_COMPILER_WRAPPER
            $mpiSkipWrapperLine = Select-String -Path $cmakeCacheFile -Pattern '^MPI_SKIP_COMPILER_WRAPPER:BOOL=' | Select-Object -First 1
            if ($mpiSkipWrapperLine) {
                $mpiSkipWrapperValue = $mpiSkipWrapperLine.Line -replace '^MPI_SKIP_COMPILER_WRAPPER:BOOL=', ''
                Write-Host "MPI_SKIP_COMPILER_WRAPPER: $mpiSkipWrapperValue" -ForegroundColor White
            } else {
                Write-Host "MPI_SKIP_COMPILER_WRAPPER: (not found in cache)" -ForegroundColor Yellow
            }
            
            # Extract MPI_Fortran_COMPILER (should be ifx.exe, not mpiifx.bat)
            $mpiFortranCompilerLine = Select-String -Path $cmakeCacheFile -Pattern '^MPI_Fortran_COMPILER:STRING=' | Select-Object -First 1
            if ($mpiFortranCompilerLine) {
                $mpiFortranCompilerValue = $mpiFortranCompilerLine.Line -replace '^MPI_Fortran_COMPILER:STRING=', ''
                Write-Host "MPI_Fortran_COMPILER: $mpiFortranCompilerValue" -ForegroundColor White
            } else {
                Write-Host "MPI_Fortran_COMPILER: (not found in cache)" -ForegroundColor Yellow
            }
            
            # Extract MPI_Fortran_LIB_NAMES (should NOT be impi - success criteria)
            $mpiFortranLibNamesLine = Select-String -Path $cmakeCacheFile -Pattern '^MPI_Fortran_LIB_NAMES:STRING=' | Select-Object -First 1
            if ($mpiFortranLibNamesLine) {
                $mpiFortranLibNamesValue = $mpiFortranLibNamesLine.Line -replace '^MPI_Fortran_LIB_NAMES:STRING=', ''
                Write-Host "MPI_Fortran_LIB_NAMES: $mpiFortranLibNamesValue" -ForegroundColor White
                if ($mpiFortranLibNamesValue -eq "impi") {
                    Write-Host "  WARNING: MPI_Fortran_LIB_NAMES is 'impi' - Intel MPI detected instead of MS-MPI!" -ForegroundColor Red
                } else {
                    Write-Host "  OK: MPI_Fortran_LIB_NAMES is not 'impi'" -ForegroundColor Green
                }
            } else {
                Write-Host "MPI_Fortran_LIB_NAMES: (not found in cache)" -ForegroundColor Yellow
            }
        }
    } else {
        Write-Host "  Warning: CMakeCache.txt not found at $cmakeCacheFile" -ForegroundColor Yellow
    }
    Write-Host ""
    
    # Workaround: Generate VendorFFTW.lib shim to avoid LNK1181
    # This is needed because the link line contains VendorFFTW.lib (not a full path),
    # even though MKL libraries are already in the link line.
    if ($FftwVendor -eq "Intel_FFTW3" -or $FftwVendor -eq "Intel_DFTI") {
        Write-Host "=== Creating VendorFFTW shim library ===" -ForegroundColor Cyan
        $shimSrcDir = Join-Path $BuildDir "vendorfftw_shim"
        $shimBuildDir = Join-Path $shimSrcDir "build"
        $targetLibPath = Join-Path $BuildDir "VendorFFTW.lib"
        
        Write-Host "Creating VendorFFTW shim lib at: $targetLibPath" -ForegroundColor Yellow
        
        # Create shim source directory
        New-Item -ItemType Directory -Force -Path $shimSrcDir | Out-Null
        New-Item -ItemType Directory -Force -Path $shimBuildDir | Out-Null
        
        # Create CMakeLists.txt for shim
        $shimCMakeLists = @"
cmake_minimum_required(VERSION 3.15)
project(VendorFFTWShim C)

file(WRITE "`${CMAKE_BINARY_DIR}/vendorfftw_dummy.c" "void vendorfftw_dummy(void){}")

add_library(VendorFFTW STATIC "`${CMAKE_BINARY_DIR}/vendorfftw_dummy.c")

set_target_properties(VendorFFTW PROPERTIES OUTPUT_NAME "VendorFFTW")
"@
        $shimCMakeListsPath = Join-Path $shimSrcDir "CMakeLists.txt"
        Set-Content -Path $shimCMakeListsPath -Value $shimCMakeLists -Encoding UTF8
        
        # Configure shim project
        $shimConfigureCmd = "$vsInitPrefix && $oneApiInit && `"$cmakePath`" -S `"$shimSrcDir`" -B `"$shimBuildDir`" -G `"$cmakeGenerator`" -DCMAKE_BUILD_TYPE=Release"
        Write-Host "  Configuring VendorFFTW shim..." -ForegroundColor Gray
        & cmd.exe /c $shimConfigureCmd
        $shimConfigureExit = $LASTEXITCODE
        if ($shimConfigureExit -ne 0) {
            throw "Failed to configure VendorFFTW shim (exit code $shimConfigureExit)"
        }
        
        # Build shim project
        $shimBuildCmd = "$vsInitPrefix && $oneApiInit && `"$cmakePath`" --build `"$shimBuildDir`" --config Release"
        Write-Host "  Building VendorFFTW shim..." -ForegroundColor Gray
        & cmd.exe /c $shimBuildCmd
        $shimBuildExit = $LASTEXITCODE
        if ($shimBuildExit -ne 0) {
            throw "Failed to build VendorFFTW shim (exit code $shimBuildExit)"
        }
        
        # Find the generated VendorFFTW.lib
        $foundLibs = Get-ChildItem -Path $shimBuildDir -Recurse -Filter "VendorFFTW.lib" -ErrorAction SilentlyContinue
        if ($foundLibs.Count -eq 0) {
            throw "VendorFFTW.lib not found in shim build directory: $shimBuildDir"
        }
        $sourceLibPath = $foundLibs[0].FullName
        Write-Host "  Found shim library at: $sourceLibPath" -ForegroundColor Green
        
        # Ensure target lib directory exists
        $targetLibDir = Split-Path -Parent $targetLibPath
        New-Item -ItemType Directory -Force -Path $targetLibDir | Out-Null
        
        # Copy to target location
        Copy-Item -Path $sourceLibPath -Destination $targetLibPath -Force
        Write-Host "  Copied to: $targetLibPath" -ForegroundColor Green
        
        # Verify the file exists
        if (-not (Test-Path $targetLibPath)) {
            throw "Failed to create VendorFFTW.lib at $targetLibPath. Shim build output: $shimBuildDir"
        }
        
        Write-Host "  VendorFFTW shim library created successfully" -ForegroundColor Green
        Write-Host ""
    }
    
    Write-Host "=== CMake Build (verbose) ===" -ForegroundColor Yellow
    
    # Capture full stdout+stderr to log file
    $buildLogFile = Join-Path $BuildDir "build.log"
    Write-Host "  Build output is being redirected to: $buildLogFile" -ForegroundColor Cyan
    Write-Host "  Check build.log for full build output" -ForegroundColor Cyan
    Write-Host ""
    
    # Use Start-Process with redirection to capture output to log file
    # This is more reliable than ProcessStartInfo with event handlers
    # Wrap the entire command in parentheses and redirect
    $buildCmdWithRedirect = "($buildCmd) > `"$buildLogFile`" 2>&1"
    
    try {
        $proc = Start-Process -FilePath "cmd.exe" -ArgumentList "/c $buildCmdWithRedirect" -NoNewWindow -Wait -PassThru
        
        $exitCode = $proc.ExitCode
    } catch {
        Write-Host "  ERROR: Failed to execute build command: $_" -ForegroundColor Red
        $exitCode = 1
        # Ensure log file exists even on error
        if (-not (Test-Path $buildLogFile)) {
            "Build command failed: $_" | Out-File -FilePath $buildLogFile -Encoding UTF8
        }
        throw
    }
    
    # Check if build succeeded
    if ($exitCode -ne 0) {
        Write-Host ""
        Write-Host "ERROR: Build failed with exit code $exitCode" -ForegroundColor Red
        Write-Host ""
        Write-Host "=== Extracting error context from build log ===" -ForegroundColor Yellow
        
        if (Test-Path $buildLogFile) {
            # Find first occurrence of "FAILED:" or "error:" with context (+-30 lines)
            $failedMatch = Select-String -Path $buildLogFile -Pattern "FAILED:" -Context 30,30 | Select-Object -First 1
            $errorMatch = Select-String -Path $buildLogFile -Pattern "error:" -Context 30,30 -CaseSensitive:$false | Select-Object -First 1
            
            if ($failedMatch) {
                Write-Host ""
                Write-Host "First FAILED: occurrence (with context):" -ForegroundColor Cyan
                Write-Host "---" -ForegroundColor Gray
                $failedMatch.Context.PreContext | ForEach-Object { Write-Host $_ }
                Write-Host $failedMatch.Line -ForegroundColor Red
                $failedMatch.Context.PostContext | ForEach-Object { Write-Host $_ }
                Write-Host "---" -ForegroundColor Gray
            } elseif ($errorMatch) {
                Write-Host ""
                Write-Host "First error: occurrence (with context):" -ForegroundColor Cyan
                Write-Host "---" -ForegroundColor Gray
                $errorMatch.Context.PreContext | ForEach-Object { Write-Host $_ }
                Write-Host $errorMatch.Line -ForegroundColor Red
                $errorMatch.Context.PostContext | ForEach-Object { Write-Host $_ }
                Write-Host "---" -ForegroundColor Gray
            } else {
                Write-Host "  No 'FAILED:' or 'error:' patterns found in log" -ForegroundColor Yellow
                Write-Host "  Showing last 60 lines of build log:" -ForegroundColor Yellow
                Write-Host "---" -ForegroundColor Gray
                Get-Content $buildLogFile -Tail 60 | ForEach-Object { Write-Host $_ }
                Write-Host "---" -ForegroundColor Gray
            }
            Write-Host ""
            Write-Host "Full build log available at: $buildLogFile" -ForegroundColor Cyan
        } else {
            Write-Host "  WARNING: Build log file not found at $buildLogFile" -ForegroundColor Yellow
        }
        
        Write-Host ""
        exit $exitCode
    } else {
        # Build succeeded - show last 60 lines
        Write-Host ""
        Write-Host "Build completed successfully. Showing last 60 lines:" -ForegroundColor Green
        Write-Host "---" -ForegroundColor Gray
        if (Test-Path $buildLogFile) {
            Get-Content $buildLogFile -Tail 60 -ErrorAction SilentlyContinue | ForEach-Object { Write-Host $_ }
        }
        Write-Host "---" -ForegroundColor Gray
        Write-Host ""
        Write-Host "Full build log available at: $buildLogFile" -ForegroundColor Cyan
    }
} finally {
    # Restore error action preference
    $ErrorActionPreference = $oldErrorAction
}

Write-Host ""
Write-Host "Step 9: Renaming executables from *.x.exe to *.exe..." -ForegroundColor Yellow
$binDir = Join-Path $BuildDir "bin"
if (Test-Path $binDir) {
    $renamedCount = 0
    Get-ChildItem -Path $binDir -Filter "*.x.exe" | ForEach-Object {
        $newName = $_.Name -replace '\.x\.exe$', '.exe'
        $newPath = Join-Path $binDir $newName
        # Always overwrite existing .exe with newly built .x.exe to prevent stale binaries
        if (Test-Path $newPath) {
            # Compare timestamps - if .x.exe is newer or same, replace
            $xExeTime = $_.LastWriteTime
            $exeTime = (Get-Item $newPath).LastWriteTime
            if ($xExeTime -ge $exeTime) {
                Remove-Item -Path $newPath -Force -ErrorAction SilentlyContinue
                Copy-Item -Path $_.FullName -Destination $newPath -Force
                Write-Host "  Replaced: $newName (with newly built $($_.Name))" -ForegroundColor Green
                $renamedCount++
            } else {
                Write-Host "  Warning: $newName is newer than $($_.Name), keeping existing" -ForegroundColor Yellow
            }
        } else {
            Copy-Item -Path $_.FullName -Destination $newPath -Force
            Write-Host "  Renamed: $($_.Name) -> $newName" -ForegroundColor Green
            $renamedCount++
        }
    }
    if ($renamedCount -gt 0) {
        Write-Host "  Processed $renamedCount executable(s)" -ForegroundColor Green
    } else {
        Write-Host "  No executables to rename (already renamed or none found)" -ForegroundColor Gray
    }
} else {
    Write-Host "  Warning: bin directory not found at $binDir" -ForegroundColor Yellow
}

Write-Host ""
Write-Host "Step 10: Verifying stack reserve in link command..." -ForegroundColor Yellow
$pwXExe = Join-Path $binDir "pw.x.exe"
if (Test-Path $pwXExe) {
    # Check if ninja build log contains the stack reserve flag
    # The verbose build output should have been shown, but we can check the actual executable
    # by looking for the flag in recent build output or checking the executable directly
    Write-Host "  Checking link command for /STACK:134217728..." -ForegroundColor Gray
    # Since we used --verbose, the link command was already shown
    # We'll verify by checking the actual executable with dumpbin (if available)
    $dumpbinPath = $null
    if ($env:VCToolsInstallDir) {
        $dumpbinPath = Join-Path $env:VCToolsInstallDir "bin\Hostx64\x64\dumpbin.exe"
        if (-not (Test-Path $dumpbinPath)) {
            $dumpbinPath = $null
        }
    }
    if (-not $dumpbinPath) {
        # Try to find dumpbin in VS paths
        $vsPaths = @(
            "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC",
            "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Tools\MSVC",
            "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Tools\MSVC"
        )
        foreach ($vsPath in $vsPaths) {
            if (Test-Path $vsPath) {
                $msvcVersions = Get-ChildItem $vsPath -Directory -ErrorAction SilentlyContinue | Sort-Object Name -Descending
                foreach ($version in $msvcVersions) {
                    $candidate = Join-Path $version.FullName "bin\Hostx64\x64\dumpbin.exe"
                    if (Test-Path $candidate) {
                        $dumpbinPath = $candidate
                        break
                    }
                }
                if ($dumpbinPath) { break }
            }
        }
    }
    
    if ($dumpbinPath) {
        try {
            $output = & $dumpbinPath /headers $pwXExe 2>&1 | Select-String -Pattern "size of stack reserve" -CaseSensitive:$false
            if ($output) {
                $line = if ($output -is [array]) { $output[0].Line } else { $output.Line }
                Write-Host "  Stack reserve: $($line.Trim())" -ForegroundColor Green
                # Parse hex value
                if ($line -match "\s*([0-9A-Fa-f]+)\s+size of stack reserve") {
                    $hexValue = $Matches[1]
                    try {
                        $bytes = [Convert]::ToInt64($hexValue, 16)
                        $mb = [math]::Round($bytes / 1MB, 2)
                        if ($mb -eq 128) {
                            Write-Host "  Verified: 128 MB stack reserve is correctly applied" -ForegroundColor Green
                        } else {
                            Write-Warning "  WARNING: Stack reserve is $mb MB, expected 128 MB"
                        }
                    } catch {
                        Write-Host "  (Could not parse hex value: $hexValue)" -ForegroundColor Yellow
                    }
                }
            } else {
                Write-Warning "  Could not find stack reserve in dumpbin output"
            }
        } catch {
            Write-Host "  (dumpbin check failed: $_)" -ForegroundColor Yellow
        }
    } else {
        Write-Host "  (dumpbin not found, skipping verification)" -ForegroundColor Gray
    }
} else {
    Write-Host "  (pw.x.exe not found, skipping verification)" -ForegroundColor Gray
}

Write-Host ""
Write-Host "==========================================" -ForegroundColor Green
Write-Host "Build completed successfully!" -ForegroundColor Green
Write-Host "==========================================" -ForegroundColor Green
Write-Host ""
Write-Host "Executables should be in: $BuildDir\bin" -ForegroundColor Cyan
Write-Host ""

exit 0

