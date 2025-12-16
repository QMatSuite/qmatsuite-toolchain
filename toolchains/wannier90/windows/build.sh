#!/bin/bash
# Wannier90 v3.1.0 build script for Windows (Intel ifort)
# This script downloads, configures, builds, and packages Wannier90
# IMPORTANT: Uses Intel Fortran Classic (ifort) from oneAPI, NOT MinGW gfortran
# MSYS2 may be used for build tools (make, tar, bash) but NOT for the compiler

set -euo pipefail

# Configuration
W90_VERSION="3.1.0"
W90_TARBALL_URL="https://github.com/wannier-developers/wannier90/archive/refs/tags/v${W90_VERSION}.tar.gz"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${BUILD_DIR:-/c/wannier90-build/wannier90-${W90_VERSION}}"
PREFIX="${PREFIX:-/c/wannier90-install/wannier90-${W90_VERSION}}"

# Print key environment variables
echo "=========================================="
echo "Wannier90 Windows Build Script (Intel ifort)"
echo "=========================================="
echo "W90_VERSION: ${W90_VERSION}"
echo "W90_TARBALL_URL: ${W90_TARBALL_URL}"
echo "SCRIPT_DIR: ${SCRIPT_DIR}"
echo "BUILD_DIR: ${BUILD_DIR}"
echo "PREFIX: ${PREFIX}"
echo "FC: ${FC:-ifort}"
echo "CC: ${CC:-icl}"
echo "=========================================="

# Verify Intel Fortran is available
if ! command -v "${FC:-ifort}" &> /dev/null; then
    echo "ERROR: Intel Fortran compiler (ifort) not found in PATH!"
    echo "Please ensure Intel oneAPI is installed and environment is set up."
    exit 1
fi

# Create build directory
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

# Download source if not already present
if [ ! -f "wannier90-${W90_VERSION}.tar.gz" ]; then
    echo "Downloading Wannier90 ${W90_VERSION} source..."
    if command -v wget &> /dev/null; then
        wget -q "${W90_TARBALL_URL}" -O "wannier90-${W90_VERSION}.tar.gz"
    elif command -v curl &> /dev/null; then
        curl -L -s "${W90_TARBALL_URL}" -o "wannier90-${W90_VERSION}.tar.gz"
    else
        echo "ERROR: Neither wget nor curl found!"
        exit 1
    fi
fi

# Extract source
if [ ! -d "wannier90-${W90_VERSION}" ]; then
    echo "Extracting source..."
    tar -xzf "wannier90-${W90_VERSION}.tar.gz"
fi

cd "wannier90-${W90_VERSION}"

# Copy make.inc from toolchain directory
echo "Copying make.inc..."
cp "${SCRIPT_DIR}/make.inc" ./make.inc

# Verify make.inc exists
if [ ! -f "./make.inc" ]; then
    echo "ERROR: make.inc not found!"
    exit 1
fi

# Show compiler info
echo "Compiler information:"
${FC:-ifort} -V || ${FC:-ifort} --version || true
which ${FC:-ifort} || true

# Build
echo "Building Wannier90..."
make clean || true
make all

# Verify build
if [ ! -f "wannier90.x" ] && [ ! -f "wannier90.exe" ]; then
    echo "ERROR: wannier90.x or wannier90.exe not built!"
    exit 1
fi

echo "Build successful! Built executables:"
ls -lh *.x *.exe 2>/dev/null || true

# Install
echo "Installing to ${PREFIX}..."
mkdir -p "${PREFIX}/bin"
mkdir -p "${PREFIX}/share/wannier90"

# Copy executables (handle both .x and .exe extensions)
cp *.x "${PREFIX}/bin/" 2>/dev/null || true
cp *.exe "${PREFIX}/bin/" 2>/dev/null || true

# Copy documentation
if [ -f "README" ]; then
    cp README "${PREFIX}/share/wannier90/"
fi
if [ -f "LICENSE" ]; then
    cp LICENSE "${PREFIX}/share/wannier90/"
fi
if [ -f "LICENCE" ]; then
    cp LICENCE "${PREFIX}/share/wannier90/LICENSE"
fi

echo "Installation complete!"
echo "Binaries: ${PREFIX}/bin/"
echo "Documentation: ${PREFIX}/share/wannier90/"
ls -R "${PREFIX}"

