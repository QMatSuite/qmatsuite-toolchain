#!/bin/bash
# Wannier90 v3.1.0 build script for macOS
# This script downloads, configures, builds, and packages Wannier90
# IMPORTANT: Uses Apple Accelerate framework with -ff2c (NOT OpenBLAS)

set -euo pipefail

# Configuration
W90_VERSION="${W90_VERSION:-3.1.0}"
W90_URL="${W90_URL:-https://github.com/wannier-developers/wannier90/archive/refs/tags/v${W90_VERSION}.tar.gz}"
W90_PREFIX="${W90_PREFIX:-${RUNNER_TEMP:-/tmp}/wannier90-install}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${BUILD_DIR:-${RUNNER_TEMP:-/tmp}/wannier90-build}"

# Print key environment variables
echo "=========================================="
echo "Wannier90 macOS Build Script"
echo "=========================================="
echo "W90_VERSION: ${W90_VERSION}"
echo "W90_URL: ${W90_URL}"
echo "W90_PREFIX: ${W90_PREFIX}"
echo "BUILD_DIR: ${BUILD_DIR}"
echo "SCRIPT_DIR: ${SCRIPT_DIR}"
echo "FC: ${FC:-gfortran}"
echo "CC: ${CC:-cc}"
echo "=========================================="

# Show compiler version
echo "Compiler information:"
${FC:-gfortran} --version || true
which ${FC:-gfortran} || true

# Create build directory
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

# Download source
TARBALL="wannier90-${W90_VERSION}.tar.gz"
if [ ! -f "${TARBALL}" ]; then
    echo "Downloading Wannier90 ${W90_VERSION} source..."
    if command -v curl &> /dev/null; then
        curl -L -s "${W90_URL}" -o "${TARBALL}"
    elif command -v wget &> /dev/null; then
        wget -q "${W90_URL}" -O "${TARBALL}"
    else
        echo "ERROR: Neither curl nor wget found!"
        exit 1
    fi
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

# Verify Accelerate framework is used (not OpenBLAS)
if grep -qi "openblas" ./make.inc; then
    echo "ERROR: make.inc contains OpenBLAS references. macOS must use Accelerate framework!"
    exit 1
fi
if ! grep -q "Accelerate" ./make.inc; then
    echo "ERROR: make.inc does not reference Accelerate framework!"
    exit 1
fi
if ! grep -q "ff2c" ./make.inc; then
    echo "WARNING: make.inc does not include -ff2c flag for Accelerate compatibility!"
fi

# Build
echo "Building Wannier90..."
make clean || true
make -j2 F90="${FC:-gfortran}"

# Verify build outputs
if [ ! -f "./wannier90.x" ]; then
    echo "ERROR: wannier90.x not built!"
    exit 1
fi

if [ ! -f "./postw90.x" ]; then
    echo "ERROR: postw90.x not built!"
    exit 1
fi

echo "Build successful! Built executables:"
ls -lh *.x || true

# Install
echo "Installing to ${W90_PREFIX}..."
mkdir -p "${W90_PREFIX}/bin"
mkdir -p "${W90_PREFIX}/share/wannier90"

# Copy executables
cp wannier90.x "${W90_PREFIX}/bin/"
cp postw90.x "${W90_PREFIX}/bin/"

# Copy any other .x executables that might exist
find . -maxdepth 1 -name "*.x" -type f -exec cp {} "${W90_PREFIX}/bin/" \;

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

