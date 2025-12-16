# Wannier90 on Windows

Windows-specific build notes and scripts for Wannier90 v3.1.0.

## Version

- **W90_VERSION**: 3.1.0
- **W90_URL**: https://github.com/wannier-developers/wannier90/archive/refs/tags/v3.1.0.tar.gz

## Dependencies

- Intel Fortran Classic (ifort) from oneAPI
- Intel MKL (included with oneAPI)
- Make (from MSYS2 or standalone)
- curl or tar (for downloading/extracting source)

**Important**: Windows builds **must** use Intel Fortran Classic (ifort) from oneAPI. MinGW gfortran is NOT used. MSYS2 may be used only for build tools like make/tar/bash, not for the compiler.

## Build Script

The `build.ps1` PowerShell script automates the build process:

```powershell
cd toolchains/wannier90/windows
.\build.ps1
```

Or with a custom prefix:
```powershell
$env:W90_PREFIX = "C:\path\to\install"
.\build.ps1
```

The script will:
1. Set up Intel Fortran compiler environment
2. Download the Wannier90 v3.1.0 source tarball
3. Extract and configure using `make.inc`
4. Build serial executables (`make -j2`)
5. Install to `$env:W90_PREFIX` (defaults to `C:\w90-install` in CI)

## Configuration

The `make.inc` file is based on `make.inc.ifort` template and configured for:
- Intel Fortran Classic compiler (ifort)
- Intel MKL (`-mkl=sequential`) for BLAS/LAPACK
- Preprocessor flag (`-fpp`) for consistency with QE Windows builds
- Serial build (MPI disabled by default)

## MPI Support

MPI support is **disabled by default**. To enable:
1. Set `COMMS=mpi` in `make.inc`
2. Set `MPIF90=mpiifort` in `make.inc`
3. Ensure Intel MPI is available

This is documented as a future extension.

## CI Strategy

Current CI strategy (see `.github/workflows/wannier90-windows-intel.yml`):

- Runner: `windows-latest`
- Install Intel Fortran Classic via `fortran-lang/setup-fortran@v1`
- Use `ifort` as the Fortran compiler
- Use Intel MKL for BLAS/LAPACK
- Build serial executables using `make -j2`
- Package artifacts with binaries in `bin/` and documentation in `share/`

## Artifact Name

- `wannier90-3.1.0-windows-intel-classic`

