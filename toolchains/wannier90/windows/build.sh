#!/bin/bash
# Wannier90 v3.1.0 build script for Windows (Intel ifort)
# Bash version - uses Intel Fortran Classic (ifort) from oneAPI
# IMPORTANT: Uses Intel Fortran Classic (ifort), NOT MinGW gfortran
# MSYS2 may be used for build tools (make, tar, bash) but NOT for the compiler

set -euo pipefail

# Configuration
W90_VERSION="${W90_VERSION:-3.1.0}"
W90_URL="${W90_URL:-https://github.com/wannier-developers/wannier90/archive/refs/tags/v${W90_VERSION}.tar.gz}"
W90_PREFIX="${W90_PREFIX:-/c/w90-install}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${BUILD_DIR:-/c/w90-build}"

# Print key environment variables
echo "=========================================="
echo "Wannier90 Windows Build Script (Intel ifort)"
echo "=========================================="
echo "W90_VERSION: ${W90_VERSION}"
echo "W90_URL: ${W90_URL}"
echo "W90_PREFIX: ${W90_PREFIX}"
echo "BUILD_DIR: ${BUILD_DIR}"
echo "SCRIPT_DIR: ${SCRIPT_DIR}"
echo "FC: ${FC:-ifort}"
echo "=========================================="

# Verify Intel Fortran is available
if ! command -v "${FC:-ifort}" &> /dev/null; then
    echo "ERROR: Intel Fortran compiler (ifort) not found in PATH!"
    echo "Please ensure Intel oneAPI is installed and environment is set up."
    exit 1
fi

# Show compiler info
echo "Compiler information:"
${FC:-ifort} -V || ${FC:-ifort} --version || true
which ${FC:-ifort} || true

# Verify it's Intel Fortran Classic (not gfortran)
VERSION_OUTPUT=$(${FC:-ifort} -V 2>&1 || true)
if ! echo "$VERSION_OUTPUT" | grep -qi "Intel.*Fortran.*Classic"; then
    echo "WARNING: Compiler may not be Intel Fortran Classic"
    echo "Version output: $VERSION_OUTPUT"
fi

# Create build directory
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

# Download source
TARBALL="wannier90-${W90_VERSION}.tar.gz"
if [ ! -f "${TARBALL}" ]; then
    echo "Downloading Wannier90 ${W90_VERSION} source..."
    curl -L -s "${W90_URL}" -o "${TARBALL}"
fi

# Extract source
SRC_DIR="${BUILD_DIR}/wannier90-${W90_VERSION}"
if [ ! -d "${SRC_DIR}" ]; then
    echo "Extracting source..."
    tar -xzf "${TARBALL}"
    
    # Find extracted directory (handle different tarball structures)
    EXTRACTED=$(find . -maxdepth 1 -type d -name "wannier90-${W90_VERSION}*" | head -1)
    if [ -z "${EXTRACTED}" ]; then
        echo "ERROR: Failed to find extracted Wannier90 source directory"
        ls -la
        exit 1
    fi
    
    # Rename to fixed name if needed
    if [ "${EXTRACTED}" != "${SRC_DIR}" ]; then
        mv "${EXTRACTED}" "${SRC_DIR}"
    fi
fi

cd "${SRC_DIR}"

# Verify Makefile exists
if [ ! -f "./Makefile" ]; then
    echo "ERROR: Makefile not found in ${SRC_DIR}"
    exit 1
fi

# Copy make.inc from toolchain directory
echo "Copying make.inc..."
cp "${SCRIPT_DIR}/make.inc" ./make.inc

# Verify make.inc exists
if [ ! -f "./make.inc" ]; then
    echo "ERROR: make.inc not found!"
    exit 1
fi

# Verify ifort is used (not gfortran)
if grep -qi "gfortran" ./make.inc; then
    echo "ERROR: make.inc contains gfortran references. Windows must use Intel ifort!"
    exit 1
fi
if ! grep -qi "ifort" ./make.inc; then
    echo "ERROR: make.inc does not reference Intel ifort!"
    exit 1
fi

# Build
echo "Building Wannier90..."
make clean || true
make -j2

# Verify build outputs
if [ ! -f "./wannier90.x" ] && [ ! -f "./wannier90.exe" ]; then
    echo "ERROR: wannier90.x or wannier90.exe not built!"
    exit 1
fi

if [ ! -f "./postw90.x" ] && [ ! -f "./postw90.exe" ]; then
    echo "ERROR: postw90.x or postw90.exe not built!"
    exit 1
fi

echo "Build successful! Built executables:"
ls -lh *.x *.exe 2>/dev/null || true

# Install
echo "Installing to ${W90_PREFIX}..."
mkdir -p "${W90_PREFIX}/bin"
mkdir -p "${W90_PREFIX}/share/wannier90"

# Copy executables (handle both .x and .exe extensions)
cp wannier90.x "${W90_PREFIX}/bin/" 2>/dev/null || true
cp wannier90.exe "${W90_PREFIX}/bin/" 2>/dev/null || true
cp postw90.x "${W90_PREFIX}/bin/" 2>/dev/null || true
cp postw90.exe "${W90_PREFIX}/bin/" 2>/dev/null || true

# Copy any other .x or .exe executables
find . -maxdepth 1 -name "*.x" -o -name "*.exe" | while read -r exe; do
    cp "$exe" "${W90_PREFIX}/bin/" 2>/dev/null || true
done

# Copy documentation
if [ -f "README.install" ]; then
    cp README.install "${W90_PREFIX}/share/wannier90/"
fi
if [ -f "README" ]; then
    cp README "${W90_PREFIX}/share/wannier90/"
fi
if [ -f "LICENSE" ]; then
    cp LICENSE "${W90_PREFIX}/share/wannier90/"
fi
if [ -f "LICENCE" ]; then
    cp LICENCE "${W90_PREFIX}/share/wannier90/LICENSE"
fi
if [ -f "COPYING" ]; then
    cp COPYING "${W90_PREFIX}/share/wannier90/"
fi

echo "Installation complete!"
echo "Installation tree:"
ls -R "${W90_PREFIX}"

