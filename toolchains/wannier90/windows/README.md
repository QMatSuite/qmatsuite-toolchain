# Wannier90 on Windows

Windows-specific build notes and scripts for Wannier90 v3.1.0.

## Version

- **W90_VERSION**: 3.1.0

## Dependencies

- Intel Fortran Classic (ifort) from oneAPI
- MSYS2 (for build tools: make, tar, bash) - **NOT for the compiler**
- Make (from MSYS2 or standalone)

**Important**: Windows builds **must** use Intel ifort classic (oneAPI). MinGW gfortran is NOT used. MSYS2 may be used only for tools like make/tar/bash, not for the compiler.

## Build Script

The `build.ps1` PowerShell script automates the build process:

```powershell
cd toolchains/wannier90/windows
.\build.ps1
```

Alternatively, if using bash (via MSYS2/Git Bash):

```bash
cd toolchains/wannier90/windows
bash build.sh
```

The script will:
1. Set up Intel Fortran compiler environment
2. Download the Wannier90 v3.1.0 source tarball
3. Extract and configure using `make.inc`
4. Build all executables
5. Install to a local prefix directory

## Configuration

The `make.inc` file is based on `make.inc.ifort` template and configured for:
- Intel Fortran Classic compiler (ifort)
- Intel MKL for BLAS/LAPACK (if available)
- Windows-specific flags and paths

## CI Strategy

Current CI strategy (see `.github/workflows/wannier90-windows-intel.yml`):

- Runner: `windows-latest`
- Install Intel Fortran Classic via `fortran-lang/setup-fortran@v1`
- Use `ifort` as the Fortran compiler
- Build all executables using `make all`
- Package artifacts with binaries in `bin/` and documentation in `share/`

## Artifact Name

- `wannier90-3.1.0-windows-intel-classic`

