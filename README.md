# QMatSuite Toolchain

This repository contains build scripts and CI workflows for the **QMatSuite**
engine/toolchain ecosystem.

It is **not** meant to be used directly by end users. Instead:

- End users should download ready-to-use binaries from:
  - a dedicated QE binaries repo (e.g. `quantum-espresso-windows-exe`, planned)
  - the QMatSuite GUI, which will download and configure the right toolchains
- This repo focuses on:
  - Building Quantum ESPRESSO (QE) on Linux, macOS, and Windows
  - Later: building related tools such as Wannier90

## Layout

- `toolchains/quantum_espresso/`
  - Platform-specific notes and scripts for QE
- `toolchains/wannier90/`
  - Reserved for Wannier90 build scripts (future)
- `docs/`
  - Additional documentation, notes, and design decisions
- `.github/workflows/`
  - GitHub Actions CI definitions for QE (and later other tools)

## Status

Currently implemented:

- Minimal CI builds of **Quantum ESPRESSO 7.5** on:
  - `ubuntu-latest`
  - `macos-14`
  - `windows-latest` (Intel oneAPI + MKL, with MS-MPI installed)

No binaries are published from this repo yet; it is intended as a backend for
other QMatSuite projects.
