# Wannier90 v3.1.0 build script for Windows (Intel ifort)
# PowerShell version - uses Intel Fortran Classic (ifort) from oneAPI
# IMPORTANT: Uses Intel Fortran Classic (ifort), NOT MinGW gfortran

$ErrorActionPreference = "Stop"

# Configuration
$W90_VERSION = if ($env:W90_VERSION) { $env:W90_VERSION } else { "3.1.0" }
$W90_URL = if ($env:W90_URL) { $env:W90_URL } else { "https://github.com/wannier-developers/wannier90/archive/refs/tags/v${W90_VERSION}.tar.gz" }
$W90_PREFIX = if ($env:W90_PREFIX) { $env:W90_PREFIX } else { "C:\w90-install" }

$SCRIPT_DIR = Split-Path -Parent $MyInvocation.MyCommand.Path
$BUILD_DIR = if ($env:BUILD_DIR) { $env:BUILD_DIR } else { "C:\w90-build" }

# Print key environment variables
Write-Host "=========================================="
Write-Host "Wannier90 Windows Build Script (Intel ifort)"
Write-Host "=========================================="
Write-Host "W90_VERSION: $W90_VERSION"
Write-Host "W90_URL: $W90_URL"
Write-Host "W90_PREFIX: $W90_PREFIX"
Write-Host "BUILD_DIR: $BUILD_DIR"
Write-Host "SCRIPT_DIR: $SCRIPT_DIR"
Write-Host "FC: $env:FC"
Write-Host "=========================================="

# Verify Intel Fortran is available
$ifort = if ($env:FC) { $env:FC } else { "ifort" }
if (-not (Get-Command $ifort -ErrorAction SilentlyContinue)) {
    Write-Host "ERROR: Intel Fortran compiler (ifort) not found in PATH!"
    Write-Host "Please ensure Intel oneAPI is installed and environment is set up."
    exit 1
}

# Show compiler info
Write-Host "Compiler information:"
& $ifort -V 2>&1 | Out-Host
Get-Command $ifort | Select-Object -ExpandProperty Source

# Verify it's Intel Fortran Classic (not gfortran)
$versionOutput = & $ifort -V 2>&1 | Out-String
if ($versionOutput -notmatch "Intel.*Fortran.*Classic") {
    Write-Host "WARNING: Compiler may not be Intel Fortran Classic"
    Write-Host "Version output: $versionOutput"
}

# Create build directory
New-Item -ItemType Directory -Force -Path $BUILD_DIR | Out-Null
Set-Location $BUILD_DIR

# Download source
$tarball = "wannier90-${W90_VERSION}.tar.gz"
if (-not (Test-Path $tarball)) {
    Write-Host "Downloading Wannier90 ${W90_VERSION} source..."
    Invoke-WebRequest -Uri $W90_URL -OutFile $tarball
}

# Extract source (use tar if available, otherwise use 7z or other)
$srcDir = Join-Path $BUILD_DIR "wannier90-${W90_VERSION}"
if (-not (Test-Path $srcDir)) {
    Write-Host "Extracting source..."
    if (Get-Command tar -ErrorAction SilentlyContinue) {
        tar -xzf $tarball
    } else {
        Write-Host "ERROR: tar command not found. Please install MSYS2 or ensure tar is in PATH."
        exit 1
    }
    
    # Find extracted directory
    $extracted = Get-ChildItem -Directory | Where-Object { $_.Name -like "wannier90-${W90_VERSION}*" } | Select-Object -First 1
    if (-not $extracted) {
        Write-Host "ERROR: Failed to find extracted Wannier90 source directory"
        Get-ChildItem
        exit 1
    }
    
    # Rename to fixed name if needed
    if ($extracted.FullName -ne $srcDir) {
        Move-Item $extracted.FullName $srcDir
    }
}

Set-Location $srcDir

# Verify Makefile exists
if (-not (Test-Path ".\Makefile")) {
    Write-Host "ERROR: Makefile not found in $srcDir"
    exit 1
}

# Copy make.inc from toolchain directory
Write-Host "Copying make.inc..."
Copy-Item "${SCRIPT_DIR}\make.inc" -Destination ".\make.inc" -Force

# Verify make.inc exists
if (-not (Test-Path ".\make.inc")) {
    Write-Host "ERROR: make.inc not found!"
    exit 1
}

# Verify ifort is used (not gfortran)
$makeIncContent = Get-Content ".\make.inc" -Raw
if ($makeIncContent -match "gfortran") {
    Write-Host "ERROR: make.inc contains gfortran references. Windows must use Intel ifort!"
    exit 1
}
if ($makeIncContent -notmatch "ifort") {
    Write-Host "ERROR: make.inc does not reference Intel ifort!"
    exit 1
}

# Build (use bash for make, as it's more reliable on Windows)
Write-Host "Building Wannier90..."
if (Get-Command bash -ErrorAction SilentlyContinue) {
    bash -c "make clean || true"
    bash -c "make -j2"
} elseif (Get-Command make -ErrorAction SilentlyContinue) {
    Write-Host "WARNING: bash not found, trying make directly..."
    make clean 2>&1 | Out-Null
    make -j2
} else {
    Write-Host "ERROR: Neither bash nor make found!"
    exit 1
}

# Verify build outputs
if ((-not (Test-Path "wannier90.x")) -and (-not (Test-Path "wannier90.exe"))) {
    Write-Host "ERROR: wannier90.x or wannier90.exe not built!"
    exit 1
}

if ((-not (Test-Path "postw90.x")) -and (-not (Test-Path "postw90.exe"))) {
    Write-Host "ERROR: postw90.x or postw90.exe not built!"
    exit 1
}

Write-Host "Build successful! Built executables:"
Get-ChildItem -Filter "*.x","*.exe" | ForEach-Object { Write-Host $_.FullName }

# Install
Write-Host "Installing to ${W90_PREFIX}..."
New-Item -ItemType Directory -Force -Path "${W90_PREFIX}\bin" | Out-Null
New-Item -ItemType Directory -Force -Path "${W90_PREFIX}\share\wannier90" | Out-Null

# Copy executables (handle both .x and .exe extensions)
Copy-Item "wannier90.x" "${W90_PREFIX}\bin\" -ErrorAction SilentlyContinue
Copy-Item "wannier90.exe" "${W90_PREFIX}\bin\" -ErrorAction SilentlyContinue
Copy-Item "postw90.x" "${W90_PREFIX}\bin\" -ErrorAction SilentlyContinue
Copy-Item "postw90.exe" "${W90_PREFIX}\bin\" -ErrorAction SilentlyContinue

# Copy any other .x or .exe executables
Get-ChildItem -Filter "*.x","*.exe" | ForEach-Object {
    Copy-Item $_.FullName "${W90_PREFIX}\bin\" -ErrorAction SilentlyContinue
}

# Copy documentation
if (Test-Path "README.install") {
    Copy-Item "README.install" "${W90_PREFIX}\share\wannier90\"
}
if (Test-Path "README") {
    Copy-Item "README" "${W90_PREFIX}\share\wannier90\"
}
if (Test-Path "LICENSE") {
    Copy-Item "LICENSE" "${W90_PREFIX}\share\wannier90\"
}
if (Test-Path "LICENCE") {
    Copy-Item "LICENCE" "${W90_PREFIX}\share\wannier90\LICENSE"
}
if (Test-Path "COPYING") {
    Copy-Item "COPYING" "${W90_PREFIX}\share\wannier90\"
}

Write-Host "Installation complete!"
Write-Host "Installation tree:"
Get-ChildItem -Recurse $W90_PREFIX | ForEach-Object { Write-Host $_.FullName }

