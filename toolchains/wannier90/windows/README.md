# Wannier90 on Windows

Windows-specific build notes and scripts for Wannier90 v3.1.0.

## Version

- **W90_VERSION**: 3.1.0
- **W90_URL**: https://github.com/wannier-developers/wannier90/archive/refs/tags/v3.1.0.tar.gz

## Current CI Strategy

Current CI strategy (see `.github/workflows/wannier90-windows-mingw.yml`):

- Runner: `windows-latest`
- Build toolchain: **MSYS2 with MinGW-w64**
- Compiler: **gfortran** (MinGW-w64)
- BLAS/LAPACK: **OpenBLAS** (MinGW-w64 package)
- Build serial executables using `make -j2`

## Artifact

The CI produces a redistributable artifact:

- **Artifact name**: `wannier90-3.1.0-windows-mingw`
- **Contents**:
  - `bin/wannier90.x.exe` - Main Wannier90 executable
  - `bin/postw90.x.exe` - Post-processing executable (if built)
  - `bin/libopenblas*.dll` - OpenBLAS runtime DLLs
  - `bin/libgfortran-*.dll` - GNU Fortran runtime DLLs
  - `bin/libquadmath-*.dll` - Quadmath library DLLs
  - `bin/libgcc_s_seh-1.dll` - GCC runtime DLL
  - `bin/libwinpthread-1.dll` - POSIX threads DLL
  - `share/wannier90/` - Documentation and license files

## Usage

The binaries are intended for standalone use and will also be integrated into QMatSuite later.

### Basic Example

```cmd
wannier90.x.exe si
```

This will run Wannier90 with the input file `si.win` (assuming standard Wannier90 naming).

### Running from Command Prompt

1. Extract the artifact to a directory (e.g., `C:\wannier90-3.1.0`)
2. Open Command Prompt
3. Navigate to the directory containing your Wannier90 input files
4. Run Wannier90:
   ```cmd
   C:\wannier90-3.1.0\bin\wannier90.x.exe si
   ```

Make sure all DLLs are in the same directory as `wannier90.x.exe`, or in a directory that is in your `PATH`.

## Configuration

The `make.inc` file is configured for:
- GNU Fortran compiler (gfortran) from MinGW-w64
- OpenBLAS for BLAS/LAPACK
- Serial build (MPI disabled by default)

## MPI Support

MPI support is **disabled by default**. To enable:
1. Set `COMMS=mpi` in `make.inc`
2. Set `MPIF90=mpif90` in `make.inc`
3. Ensure MPI libraries are available in MSYS2

This is documented as a future extension.

## Build Requirements

To build Wannier90 on Windows locally using the same toolchain:

1. Install MSYS2 from https://www.msys2.org/
2. Open MSYS2 MinGW 64-bit terminal
3. Install required packages:
   ```bash
   pacman -S mingw-w64-x86_64-gcc
   pacman -S mingw-w64-x86_64-gcc-fortran
   pacman -S mingw-w64-x86_64-openblas
   pacman -S mingw-w64-x86_64-fftw
   pacman -S make tar wget
   ```
4. Follow the CI workflow steps in `.github/workflows/wannier90-windows-mingw.yml`

## Patches

A patch from the [MSYS2 MINGW-packages repository](https://github.com/msys2/MINGW-packages/tree/master/mingw-w64-wannier90) is applied to fix executable suffix handling for Windows/MinGW builds. This patch is stored in `toolchains/wannier90/windows/patches/` and is automatically applied during the CI build process.

The patch makes the executable suffix configurable via an `EXE_EXT` variable (defaults to `.x`), allowing flexibility in how executables are named on Windows.

See `toolchains/wannier90/windows/patches/README.md` for details about the patch.

## Migration from Intel oneAPI

This toolchain replaces the previous Intel oneAPI/ifort based Windows builds. The MinGW-w64 toolchain provides:

- Open-source compiler (gfortran)
- No proprietary dependencies
- Better compatibility with open-source scientific software
- Easier redistribution (no Intel license restrictions)
