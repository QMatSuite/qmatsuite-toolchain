# Wannier90 on Linux

Linux-specific build notes and scripts for Wannier90 v3.1.0.

## Version

- **W90_VERSION**: 3.1.0
- **W90_URL**: https://github.com/wannier-developers/wannier90/archive/refs/tags/v3.1.0.tar.gz

## Dependencies

- GNU Fortran compiler (gfortran)
- BLAS/LAPACK libraries (system-provided)
- Make
- curl or wget (for downloading source)

Install on Ubuntu/Debian:
```bash
sudo apt-get update
sudo apt-get install -y build-essential gfortran make wget curl \
  libopenblas-dev liblapack-dev
```

## Build Script

The `build.sh` script automates the build process:

```bash
cd toolchains/wannier90/linux
./build.sh
```

Or with a custom prefix:
```bash
W90_PREFIX=/path/to/install ./build.sh
```

The script will:
1. Download the Wannier90 v3.1.0 source tarball
2. Extract and configure using `make.inc`
3. Build serial executables (`make -j2`)
4. Install to `$W90_PREFIX` (defaults to `$RUNNER_TEMP/wannier90-install` in CI)

## Configuration

The `make.inc` file is based on `make.inc.gfort` template and configured for:
- GNU Fortran compiler (gfortran)
- System BLAS/LAPACK (`-llapack -lblas`)
- Serial build (MPI disabled by default)
- Standard optimization flags (`-O3`)

## MPI Support

MPI support for parallel postw90 is **disabled by default**. To enable:
1. Set `COMMS=mpi` in `make.inc`
2. Set `MPIF90=mpif90` in `make.inc`
3. Ensure MPI libraries are available

This is documented as a future extension.

## CI Strategy

Current CI strategy (see `.github/workflows/wannier90-ubuntu.yml`):

- Runner: `ubuntu-latest`
- Dependencies installed via `apt`
- Wannier90 v3.1.0 source downloaded from GitHub releases
- Configure with system BLAS/LAPACK
- Build serial executables using `make -j2`
- Package artifacts with binaries in `bin/` and documentation in `share/`

## Artifact Name

- `wannier90-3.1.0-ubuntu-latest`

