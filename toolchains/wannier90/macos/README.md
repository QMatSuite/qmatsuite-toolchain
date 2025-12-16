# Wannier90 on macOS

macOS-specific build notes and scripts for Wannier90 v3.1.0.

## Version

- **W90_VERSION**: 3.1.0
- **W90_URL**: https://github.com/wannier-developers/wannier90/archive/refs/tags/v3.1.0.tar.gz

## Dependencies

- GNU Fortran compiler (gfortran) from Homebrew
- Apple Accelerate framework (system BLAS/LAPACK)
- Make
- curl or wget (for downloading source)

Install via Homebrew:
```bash
brew update
brew install make wget curl
```

Install Fortran compiler:
```bash
# Use fortran-lang/setup-fortran action in CI, or manually:
brew install gcc
```

**Important**: macOS builds **must** use Apple Accelerate framework (`-framework Accelerate`) for BLAS/LAPACK. The build uses `-ff2c` flag for Accelerate compatibility, as recommended by the macports template.

## Build Script

The `build.sh` script automates the build process:

```bash
cd toolchains/wannier90/macos
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

The `make.inc` file is based on `make.inc.macosx.macports` template and configured for:
- GNU Fortran compiler (gfortran)
- **Apple Accelerate framework** for BLAS/LAPACK (not OpenBLAS)
- **`-ff2c` flag** for Accelerate compatibility
- Serial build (MPI disabled by default)

## MPI Support

MPI support is **disabled by default** on macOS. The upstream template notes that `-ff2c` can conflict with OpenMPI wrappers. To enable MPI:
1. Remove `-ff2c` from `FCOPTS` and `LDOPTS`
2. Set `COMMS=mpi` in `make.inc`
3. Set `MPIF90=mpif90` in `make.inc`
4. Ensure MPI libraries are available

This is documented as a future extension.

## CI Strategy

Current CI strategy (see `.github/workflows/wannier90-macos.yml`):

- Runner: `macos-14`
- Dependencies installed via Homebrew
- Fortran compiler setup via `fortran-lang/setup-fortran@v1`
- BLAS/LAPACK from Apple's **Accelerate** framework (with `-ff2c`)
- Build serial executables using `make -j2`
- Package artifacts with binaries in `bin/` and documentation in `share/`

## Artifact Name

- `wannier90-3.1.0-macos-14`

