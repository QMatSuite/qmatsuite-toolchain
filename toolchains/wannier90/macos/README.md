# Wannier90 on macOS

macOS-specific build notes and scripts for Wannier90 v3.1.0.

## Version

- **W90_VERSION**: 3.1.0

## Dependencies

- GNU Fortran compiler (gfortran) from Homebrew
- Apple Accelerate framework (system BLAS/LAPACK)
- Make
- wget or curl (for downloading source)

Install via Homebrew:
```bash
brew update
brew install gcc gfortran wget make
```

**Important**: macOS builds **must** use Apple Accelerate framework (`-framework Accelerate`) for BLAS/LAPACK. OpenBLAS is NOT used as the default.

## Build Script

The `build.sh` script automates the build process:

```bash
cd toolchains/wannier90/macos
./build.sh
```

The script will:
1. Download the Wannier90 v3.1.0 source tarball
2. Extract and configure using `make.inc`
3. Build all executables
4. Install to a local prefix directory

## Configuration

The `make.inc` file is based on `make.inc.macosx.homebrew` template and configured for:
- GNU Fortran compiler from Homebrew
- **Apple Accelerate framework** for BLAS/LAPACK (not OpenBLAS)
- Standard optimization flags

## CI Strategy

Current CI strategy (see `.github/workflows/wannier90-macos.yml`):

- Runner: `macos-14`
- Dependencies installed via Homebrew
- BLAS/LAPACK from Apple's **Accelerate** framework
- Build all executables using `make all`
- Package artifacts with binaries in `bin/` and documentation in `share/`

## Artifact Name

- `wannier90-3.1.0-macos-14`

