# Wannier90 Toolchain

This directory contains platform-specific build scripts and configuration files for
**Wannier90 v3.1.0** used in the QMatSuite ecosystem.

## Version

- **W90_VERSION**: 3.1.0
- **Source tarball**: https://github.com/wannier-developers/wannier90/archive/refs/tags/v3.1.0.tar.gz

## Subdirectories

- `linux/`   – build scripts and configuration for Linux (Ubuntu)
- `macos/`   – build scripts and configuration for macOS
- `windows/` – build scripts and configuration for Windows (Intel ifort)

## Artifact Layout

The built binaries follow a QE-like layout:

```
prefix/
  bin/
    wannier90.x
    postw90.x
    w90chk2chk.x
    w90vdw.x
    w90spn2spn.x
    w90unk2unk.x
    w90pov.x
  share/
    wannier90/
      README
      LICENSE
```

## Platform-Specific Notes

### Linux
- Uses OpenBLAS for BLAS/LAPACK
- Built with GNU Fortran (gfortran)
- See `linux/README.md` for details

### macOS
- **Must use Apple Accelerate framework** (`-framework Accelerate`) for BLAS/LAPACK
- OpenBLAS is NOT used as the default on macOS
- Built with GNU Fortran (gfortran) from Homebrew
- See `macos/README.md` for details

### Windows
- Uses Intel Fortran Classic (ifort) from oneAPI
- MSYS2 may be used for build tools (make, tar, bash) but NOT for the compiler
- MinGW gfortran is NOT used
- See `windows/README.md` for details

## CI Workflows

GitHub Actions workflows are located in `.github/workflows/`:
- `wannier90-ubuntu.yml` – Linux builds
- `wannier90-macos.yml` – macOS builds
- `wannier90-windows-intel.yml` – Windows builds with Intel ifort

Each workflow is path-scoped to trigger only when its corresponding toolchain directory changes.
