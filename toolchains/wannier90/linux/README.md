# Wannier90 on Linux

Linux-specific build notes and scripts for Wannier90 v3.1.0.

## Version

- **W90_VERSION**: 3.1.0

## Dependencies

- GNU Fortran compiler (gfortran)
- OpenBLAS (libopenblas-dev)
- LAPACK (liblapack-dev)
- Make
- wget or curl (for downloading source)

Install on Ubuntu/Debian:
```bash
sudo apt-get update
sudo apt-get install -y build-essential gfortran libopenblas-dev liblapack-dev wget make
```

## Build Script

The `build.sh` script automates the build process:

```bash
cd toolchains/wannier90/linux
./build.sh
```

The script will:
1. Download the Wannier90 v3.1.0 source tarball
2. Extract and configure using `make.inc`
3. Build all executables
4. Install to a local prefix directory

## Configuration

The `make.inc` file is based on `make.inc.gfort` template and configured for:
- GNU Fortran compiler
- OpenBLAS for BLAS/LAPACK
- Standard optimization flags

## CI Strategy

Current CI strategy (see `.github/workflows/wannier90-ubuntu.yml`):

- Runner: `ubuntu-latest`
- Dependencies installed via `apt`
- Wannier90 v3.1.0 source downloaded from GitHub releases
- Configure with OpenBLAS
- Build all executables using `make all`
- Package artifacts with binaries in `bin/` and documentation in `share/`

## Artifact Name

- `wannier90-3.1.0-ubuntu-latest`

