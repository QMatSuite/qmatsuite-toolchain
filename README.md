# QMatSuite Toolchain

This repository contains build scripts and CI workflows for the **QMatSuite**
engine/toolchain ecosystem.

Everyone is welcome to use these compilation recipes. This repository is both (1) the build infrastructure behind QMatSuite, and (2) a curated, battle-tested knowledge base for compiling Quantum ESPRESSO and Wannier90 across platforms.

It is especially valuable for Windows, where building scientific Fortran/C toolchains is significantly harder than on Linux/macOS due to compiler/runtime differences, library availability, ABI and packaging quirks, stack/heap behavior, and multi-threading edge cases. The upside is that if you reproduce the recipes here, you can obtain high-performance, native executables comparable to "best-case" builds on other platforms.

**For ready-to-run executables:** See the appropriate releases or downloads location. **For learning, reproducibility, auditing, or customization:** This repository is the canonical source.

## What This Repository Provides

This repository documents real-world compilation journeys across Linux, macOS, and Windows (MinGW, Intel oneAPI, MS-MPI). The platform-specific challenges encountered are intentional learning points, not hidden failures. Each subdirectory contains detailed READMEs documenting the complete build process, troubleshooting, and toolchain-specific considerations.

**Quantum ESPRESSO 7.5** builds are supported on:
- Linux (Ubuntu)
- macOS
- Windows (MinGW gfortran, Intel oneAPI ifort/ifx)

**Wannier90 v3.1.0** builds are supported on:
- Linux (Ubuntu)
- macOS
- Windows (MinGW, Intel oneAPI)

**Note:** Wannier90 on Windows, especially with Intel oneAPI, required CMake and development branch features not present in the stable v3.1.0 release. This is documented in the platform-specific READMEs.

## Repository Structure

- `toolchains/quantum_espresso/`
  - Platform-specific build guides, patches, and scripts for QE
- `toolchains/wannier90/`
  - Platform-specific build guides and scripts for Wannier90
- `docs/`
  - Additional documentation, notes, and design decisions
- `.github/workflows/`
  - GitHub Actions CI definitions for automated builds

CI workflows exist for continuous integration, but local builds are equally supported and documented. The CI serves as both validation and a reference implementation.

## Build Status

Currently implemented and validated:

- **Quantum ESPRESSO 7.5** builds on:
  - `ubuntu-latest`
  - `macos-14`
  - `windows-latest` (Intel oneAPI + MKL, with MS-MPI installed)

- **Wannier90 v3.1.0** builds on:
  - `ubuntu-latest` (serial build, system BLAS/LAPACK)
  - `macos-14` (serial build, Apple Accelerate framework with -ff2c)
  - `windows-latest` (serial build, Intel ifort classic + MKL)

This repository serves as both the build infrastructure for QMatSuite projects and a standalone reference for reproducible builds across platforms.

## Windows Platform Complexity

**Important:** Windows is not a single platform in this context. Builds are split by toolchain, each with distinct characteristics and requirements.

### Toolchain Split

- **MinGW (gfortran)**: Open-source compiler, requires additional patches and careful configuration
- **Intel oneAPI (ifort/ifx)**: Proprietary compiler, generally more consistent but still non-trivial

### Why Windows Is Harder

Windows compilation presents unique challenges compared to Linux/macOS:

- **Compiler/runtime differences**: MinGW gfortran and Intel oneAPI handle Fortran features differently, requiring toolchain-specific fixes
- **Library availability**: Different packaging systems and library locations
- **ABI and packaging quirks**: DLL dependencies, stack/heap behavior, and runtime linking
- **Multi-threading edge cases**: OpenMP behavior differs, requiring careful flag selection

These challenges are documented in detail in the platform-specific READMEs. When addressed correctly, Windows builds achieve high performance and native execution comparable to Linux/macOS builds.

### Platform-Specific Documentation

Detailed build instructions, troubleshooting, patches, compiler flags, and toolchain-specific notes are in:

- **Quantum ESPRESSO:**
  - [`toolchains/quantum_espresso/windows/mingw/README.md`](toolchains/quantum_espresso/windows/mingw/README.md) — MinGW (gfortran) build guide
  - [`toolchains/quantum_espresso/windows/oneapi/README.md`](toolchains/quantum_espresso/windows/oneapi/README.md) — Intel oneAPI build guide
  - [`toolchains/quantum_espresso/linux/README.md`](toolchains/quantum_espresso/linux/README.md) — Linux build guide
  - [`toolchains/quantum_espresso/macos/README.md`](toolchains/quantum_espresso/macos/README.md) — macOS build guide

- **Wannier90:**
  - [`toolchains/wannier90/windows/mingw/README.md`](toolchains/wannier90/windows/mingw/README.md) — MinGW build guide
  - [`toolchains/wannier90/windows/oneapi/README.md`](toolchains/wannier90/windows/oneapi/README.md) — Intel oneAPI build guide
  - [`toolchains/wannier90/linux/README.md`](toolchains/wannier90/linux/README.md) — Linux build guide
  - [`toolchains/wannier90/macos/README.md`](toolchains/wannier90/macos/README.md) — macOS build guide

Each README contains the complete journey from initial setup through troubleshooting, including all patches, compiler flags, and build commands needed for that platform/toolchain combination.