# QE on Windows

Windows-specific build notes and scripts for Quantum ESPRESSO.

## Current CI Strategy

Current CI strategy (see `.github/workflows/qe-windows-mingw.yml`):

- Runner: `windows-latest`
- Build toolchain: **MSYS2 with MinGW-w64**
- Compiler: **gfortran** (MinGW-w64)
- BLAS/LAPACK: **OpenBLAS** (MinGW-w64 package)
- FFTW: **FFTW3** (MinGW-w64 package)
- Build at least `pw.x` (serial build)

## Artifact

The CI produces a redistributable artifact:

- **Artifact name**: `qe-qe-7.5-windows-mingw`
- **Contents**:
  - `bin/pw.x.exe` - Main QE executable
  - `bin/libopenblas*.dll` - OpenBLAS runtime DLLs
  - `bin/libgfortran-*.dll` - GNU Fortran runtime DLLs
  - `bin/libquadmath-*.dll` - Quadmath library DLLs
  - `bin/libgcc_s_seh-1.dll` - GCC runtime DLL
  - `bin/libwinpthread-1.dll` - POSIX threads DLL

## Usage

The binaries are intended for standalone use and will also be integrated into QMatSuite later.

### Basic Example

```cmd
pw.x.exe < si_scf.in > si_scf.out
```

Make sure all DLLs are in the same directory as `pw.x.exe`, or in a directory that is in your `PATH`.

### Running from Command Prompt

1. Extract the artifact to a directory (e.g., `C:\qe-7.5`)
2. Open Command Prompt
3. Navigate to the `bin` directory:
   ```cmd
   cd C:\qe-7.5\bin
   ```
4. Run QE:
   ```cmd
   pw.x.exe < input.in > output.out
   ```

## Build Requirements

To build QE on Windows locally using the same toolchain:

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
4. Follow the CI workflow steps in `.github/workflows/qe-windows-mingw.yml`

## Migration from Intel oneAPI

This toolchain replaces the previous Intel oneAPI/ifort/MS-MPI based Windows builds. The MinGW-w64 toolchain provides:

- Open-source compiler (gfortran)
- No proprietary dependencies
- Better compatibility with open-source scientific software
- Easier redistribution (no Intel license restrictions)
