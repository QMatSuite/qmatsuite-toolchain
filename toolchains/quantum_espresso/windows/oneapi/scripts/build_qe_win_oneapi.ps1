# QE Windows oneAPI build script
param(
    [switch]$NoMpi,
    [string]$QeSourceDir = $null
)

$ErrorActionPreference = 'Stop'
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$OneApiRoot = Split-Path $ScriptDir -Parent

if (-not $QeSourceDir) {
    $QeSourceDir = Join-Path $OneApiRoot "upstream/qe"
}
if ([System.IO.Path]::IsPathRooted($QeSourceDir)) {
    $QESource = Resolve-Path $QeSourceDir -ErrorAction Stop
} else {
    $QESource = Resolve-Path (Join-Path $OneApiRoot $QeSourceDir) -ErrorAction Stop
}
$BuildDir = Join-Path $QESource "build-win-oneapi"

Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "Quantum ESPRESSO Windows oneAPI Build" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "QE source : $QESource"
Write-Host "Build dir : $BuildDir"
Write-Host ""

if (-not (Test-Path $QESource)) { throw "QE source not found: $QESource" }
if (-not (Test-Path (Join-Path $QESource "CMakeLists.txt"))) { throw "Not a QE source (CMakeLists.txt missing): $QESource" }

$requiredSubmodules = @(
    @{Path = "external/wannier90"; CheckFile = "src/parameters.F90"},
    @{Path = "external/mbd"; CheckFile = "CMakeLists.txt"},
    @{Path = "external/devxlib"; CheckFile = "src/deviceXlib_mod.f90"}
)
foreach ($sub in $requiredSubmodules) {
    $check = Join-Path (Join-Path $QESource $sub.Path) $sub.CheckFile
    if (-not (Test-Path $check)) {
        throw "Missing submodule content: $($sub.Path). Run refresh_qe_source.ps1 first."
    }
}

# oneAPI discovery
if ($env:ONEAPI_ROOT) {
    $oneApiRoot = $env:ONEAPI_ROOT
    Write-Host "Using ONEAPI_ROOT from environment: $oneApiRoot"
} else {
    $oneApiRoot = "C:\Program Files (x86)\Intel\oneAPI"
    Write-Host "Trying default ONEAPI_ROOT: $oneApiRoot"
}
$setvarsBat = Join-Path $oneApiRoot "setvars.bat"
if (-not (Test-Path $setvarsBat)) { throw "setvars.bat not found at $setvarsBat" }

# VS detection
$vsCandidate = "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\Tools\VsDevCmd.bat"
if (-not (Test-Path $vsCandidate)) { throw "VsDevCmd.bat not found; install VS Build Tools 2022." }

# CMake/Ninja
$cmake = Get-Command cmake -ErrorAction SilentlyContinue
if (-not $cmake) { $cmake = "C:\Program Files\CMake\bin\cmake.exe" }
if (-not (Test-Path $cmake)) { throw "cmake not found; install CMake." }

$ninja = Get-Command ninja -ErrorAction SilentlyContinue
if (-not $ninja) { $ninja = "ninja" } # rely on PATH in setvars/VS

# Clean build dir cache
if (Test-Path (Join-Path $BuildDir "CMakeCache.txt")) { Remove-Item (Join-Path $BuildDir "CMakeCache.txt") -Force }
if (Test-Path (Join-Path $BuildDir "CMakeFiles")) { Remove-Item (Join-Path $BuildDir "CMakeFiles") -Recurse -Force }
New-Item -ItemType Directory -Force -Path $BuildDir | Out-Null

$mpiFlag = if ($NoMpi) { "OFF" } else { "ON" }

$cmakeConfigure = @(
    "`"$cmake`"", "-S", "`"$QESource`"", "-B", "`"$BuildDir`"",
    "-G", "Ninja",
    "-DCMAKE_BUILD_TYPE=Release",
    "-DCMAKE_Fortran_COMPILER=`"$($oneApiRoot)\\compiler\\latest\\bin\\ifx.exe`"",
    "-DCMAKE_C_COMPILER=`"$($oneApiRoot)\\compiler\\latest\\bin\\icx.exe`"",
    "-DCMAKE_Fortran_COMPILER_FORCED=TRUE",
    "-DQE_CPP=`"$($oneApiRoot)\\compiler\\latest\\bin\\icx.exe`"",
    "-DQE_ENABLE_MPI=$mpiFlag",
    "-DQE_LAPACK_INTERNAL=OFF",
    "-DQE_ENABLE_OPENMP=OFF",
    "-DQE_ENABLE_SCALAPACK=OFF",
    "-DQE_ENABLE_HDF5=OFF",
    "-DQE_ENABLE_LIBXC=OFF",
    "-DQE_ENABLE_ELPA=OFF",
    "-DQE_ENABLE_TEST=OFF",
    "-DQE_WANNIER90_INTERNAL=ON",
    "-DQE_MBD_INTERNAL=ON",
    "-DQE_DEVICEXLIB_INTERNAL=ON",
    "-DQE_ENABLE_ENVIRON=NO",
    "-DQE_FFTW_VENDOR=Internal"
) -join " "

$cmakeBuild = @(
    "`"$cmake`"", "--build", "`"$BuildDir`"", "--target", "pw", "--config", "Release", "--parallel"
) -join " "

$cmdLine = "`"$vsCandidate`" -arch=x64 -host_arch=x64 >nul 2>&1 && `"$setvarsBat`" intel64 >nul 2>&1 && $cmakeConfigure && $cmakeBuild"

Write-Host "Command to execute:"
Write-Host $cmdLine

& cmd /c $cmdLine
if ($LASTEXITCODE -ne 0) { throw "Build failed with exit code $LASTEXITCODE" }

Write-Host ""
Write-Host "=========================================="
Write-Host "Build completed successfully!"
Write-Host "Executables should be in: $BuildDir\\bin"
Write-Host "=========================================="

