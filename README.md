# QMatSuite Toolchain

This repository contains build scripts and CI workflows for the **QMatSuite**
engine/toolchain ecosystem.

It is **not** meant to be used directly by end users. Instead:

- End users should download ready-to-use binaries from:
  - a dedicated QE binaries repo (e.g. `quantum-espresso-windows-exe`, planned)
  - the QMatSuite GUI, which will download and configure the right toolchains
- This repo focuses on:
  - Building Quantum ESPRESSO (QE) on Linux, macOS, and Windows
  - Building Wannier90 v3.1.0 on Linux, macOS, and Windows

## Layout

- `toolchains/quantum_espresso/`
  - Platform-specific notes and scripts for QE
- `toolchains/wannier90/`
  - Platform-specific notes and scripts for Wannier90 v3.1.0
- `docs/`
  - Additional documentation, notes, and design decisions
- `.github/workflows/`
  - GitHub Actions CI definitions for QE and Wannier90

## Status

Currently implemented:

- Minimal CI builds of **Quantum ESPRESSO 7.5** on:
  - `ubuntu-latest`
  - `macos-14`
  - `windows-latest` (Intel oneAPI + MKL, with MS-MPI installed)

- Minimal CI builds of **Wannier90 v3.1.0** on:
  - `ubuntu-latest` (serial build, system BLAS/LAPACK)
  - `macos-14` (serial build, Apple Accelerate framework with -ff2c)
  - `windows-latest` (serial build, Intel ifort classic + MKL)

No binaries are published from this repo yet; it is intended as a backend for
other QMatSuite projects.
