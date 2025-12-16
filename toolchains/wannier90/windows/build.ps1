# Wannier90 v3.1.0 build script for Windows (Intel ifort)
# PowerShell version - uses Intel Fortran Classic (ifort) from oneAPI
# IMPORTANT: Uses Intel Fortran Classic (ifort), NOT MinGW gfortran

$ErrorActionPreference = "Stop"

# Configuration
$W90_VERSION = "3.1.0"
$W90_TARBALL_URL = "https://github.com/wannier-developers/wannier90/archive/refs/tags/v${W90_VERSION}.tar.gz"
$SCRIPT_DIR = Split-Path -Parent $MyInvocation.MyCommand.Path
$BUILD_DIR = if ($env:BUILD_DIR) { $env:BUILD_DIR } else { "C:\wannier90-build\wannier90-${W90_VERSION}" }
$PREFIX = if ($env:PREFIX) { $env:PREFIX } else { "C:\wannier90-install\wannier90-${W90_VERSION}" }

# Print key environment variables
Write-Host "=========================================="
Write-Host "Wannier90 Windows Build Script (Intel ifort)"
Write-Host "=========================================="
Write-Host "W90_VERSION: $W90_VERSION"
Write-Host "W90_TARBALL_URL: $W90_TARBALL_URL"
Write-Host "SCRIPT_DIR: $SCRIPT_DIR"
Write-Host "BUILD_DIR: $BUILD_DIR"
Write-Host "PREFIX: $PREFIX"
Write-Host "FC: $env:FC"
Write-Host "=========================================="

# Verify Intel Fortran is available
$ifort = if ($env:FC) { $env:FC } else { "ifort" }
if (-not (Get-Command $ifort -ErrorAction SilentlyContinue)) {
    Write-Host "ERROR: Intel Fortran compiler (ifort) not found in PATH!"
    Write-Host "Please ensure Intel oneAPI is installed and environment is set up."
    exit 1
}

# Create build directory
New-Item -ItemType Directory -Force -Path $BUILD_DIR | Out-Null
Set-Location $BUILD_DIR

# Download source if not already present
$tarball = "wannier90-${W90_VERSION}.tar.gz"
if (-not (Test-Path $tarball)) {
    Write-Host "Downloading Wannier90 ${W90_VERSION} source..."
    Invoke-WebRequest -Uri $W90_TARBALL_URL -OutFile $tarball
}

# Extract source
$extractedDir = "wannier90-${W90_VERSION}"
if (-not (Test-Path $extractedDir)) {
    Write-Host "Extracting source..."
    tar -xzf $tarball
}

Set-Location $extractedDir

# Copy make.inc from toolchain directory
Write-Host "Copying make.inc..."
Copy-Item "${SCRIPT_DIR}\make.inc" -Destination ".\make.inc" -Force

# Verify make.inc exists
if (-not (Test-Path ".\make.inc")) {
    Write-Host "ERROR: make.inc not found!"
    exit 1
}

# Show compiler info
Write-Host "Compiler information:"
& $ifort -V 2>&1 | Out-Host
Get-Command $ifort | Select-Object -ExpandProperty Source

# Build (use bash for make, as it's more reliable on Windows)
Write-Host "Building Wannier90..."
if (Get-Command bash -ErrorAction SilentlyContinue) {
    bash -c "make clean || true"
    bash -c "make all"
} else {
    Write-Host "WARNING: bash not found, trying make directly..."
    make clean 2>&1 | Out-Null
    make all
}

# Verify build
if ((-not (Test-Path "wannier90.x")) -and (-not (Test-Path "wannier90.exe"))) {
    Write-Host "ERROR: wannier90.x or wannier90.exe not built!"
    exit 1
}

Write-Host "Build successful! Built executables:"
Get-ChildItem -Filter "*.x","*.exe" | ForEach-Object { Write-Host $_.FullName }

# Install
Write-Host "Installing to ${PREFIX}..."
New-Item -ItemType Directory -Force -Path "${PREFIX}\bin" | Out-Null
New-Item -ItemType Directory -Force -Path "${PREFIX}\share\wannier90" | Out-Null

# Copy executables
Copy-Item "*.x" "${PREFIX}\bin\" -ErrorAction SilentlyContinue
Copy-Item "*.exe" "${PREFIX}\bin\" -ErrorAction SilentlyContinue

# Copy documentation
if (Test-Path "README") {
    Copy-Item "README" "${PREFIX}\share\wannier90\"
}
if (Test-Path "LICENSE") {
    Copy-Item "LICENSE" "${PREFIX}\share\wannier90\"
}
if (Test-Path "LICENCE") {
    Copy-Item "LICENCE" "${PREFIX}\share\wannier90\LICENSE"
}

Write-Host "Installation complete!"
Write-Host "Binaries: ${PREFIX}\bin\"
Write-Host "Documentation: ${PREFIX}\share\wannier90\"
Get-ChildItem -Recurse $PREFIX | ForEach-Object { Write-Host $_.FullName }

