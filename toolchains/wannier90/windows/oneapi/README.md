# Wannier90 Windows Build (Intel oneAPI + MKL)

This directory contains scripts to build Wannier90 on Windows using Intel oneAPI compilers (ifx) and Intel MKL via CMake.

## Structure

```
toolchains/wannier90/windows/oneapi/
├── scripts/
│   ├── refresh_wannier90_source.ps1  # Clone/update Wannier90 from git
│   └── build_wannier90_win_oneapi.ps1 # Build Wannier90 with ifx + MKL (CMake)
└── upstream/                          # Wannier90 source (cloned from git)
    └── wannier90/
```

## Prerequisites

1. **Intel oneAPI HPC Toolkit** (2024 or later)
   - Includes `ifx` (Fortran) compiler
   - Includes Intel MKL (BLAS/LAPACK)
   - Download from: https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html

2. **Visual Studio 2022 Build Tools** (or full Visual Studio)
   - Required for the MSVC linker (Intel compilers need it on Windows)
   - Install "Desktop development with C++" workload

3. **CMake** (3.25+) and **Ninja**
   - CMake: Download from https://cmake.org/download/ or `choco install cmake`
   - Ninja: `choco install ninja` or download from https://github.com/ninja-build/ninja/releases

4. **Git** (for cloning Wannier90)
   - Download from: https://git-scm.com/download/win
   - Or install via Chocolatey: `choco install git`

## Quick Start

### Step 1: Clone Wannier90 source

The refresh script defaults to the `develop` branch which includes CMake support:

```powershell
.\scripts\refresh_wannier90_source.ps1
```

This clones Wannier90 from GitHub into `upstream/wannier90`.

**Note:** The script defaults to `develop` branch for CMake support. If you need a specific version (e.g., `v3.1.0` for Make build), use:
```powershell
.\scripts\refresh_wannier90_source.ps1 -W90Version v3.1.0
```

### Step 2: Build Wannier90

```powershell
.\scripts\build_wannier90_win_oneapi.ps1
```

The build script will:
1. Verify Wannier90 source has CMakeLists.txt (develop branch)
2. Discover Intel oneAPI installation
3. Set up MKL environment for CMake
4. Check for Visual Studio (for linker)
5. Check for CMake and Ninja
6. Configure and build Wannier90 using CMake
7. Output executables in `build-win-oneapi/` directory

**Note:** This build uses CMake and follows the same approach as Quantum ESPRESSO. It requires the `develop` branch which includes CMake support.

## Where Does CMakeLists.txt Come From?

The `CMakeLists.txt` file is part of Wannier90's source code in the `develop` branch on GitHub:
- **Repository**: https://github.com/wannier-developers/wannier90
- **Branch**: `develop` (contains CMake support)
- **Location**: Root of the repository (`CMakeLists.txt`)

The refresh script clones this repository and checks out the `develop` branch by default, which includes the CMake build system. The `v3.1.0` release tag does not include CMake support (only Make).

## Why CMake Instead of Make?

This build system uses CMake instead of the traditional Make-based build for several important reasons related to Windows and Intel oneAPI compatibility:

### Issues Encountered with Make Build

1. **Linker Flag Conflicts**
   - The Make build system passes compiler flags (like `-warn`) directly to the linker
   - On Windows, Intel's `ifx` compiler uses the MSVC linker (`link.exe`), which doesn't understand Fortran compiler-specific flags
   - This caused errors like: `link: unknown option -- s` when `-warn` was passed to the linker
   - **Solution**: CMake properly separates compiler flags from linker flags, avoiding this issue

2. **Flag Separation Problems**
   - Make's `make.inc` typically sets `LDOPTS=$(FCOPTS)`, copying all compiler flags to the linker
   - Flags like `-warn`, `-fpp` (preprocessor), and other compiler-only flags are not valid for the linker
   - Manual separation of flags in `make.inc` is error-prone and requires Windows-specific knowledge
   - **Solution**: CMake automatically handles flag separation based on compiler and platform

3. **Windows Path and Environment Handling**
   - Make builds require careful handling of Windows paths (backslashes vs forward slashes)
   - Environment variable propagation between PowerShell, cmd.exe, and make can be problematic
   - Unix-style commands in Makefiles (like `rm -f`) don't work natively on Windows
   - **Solution**: CMake abstracts platform differences and handles paths correctly on all platforms

4. **MKL Library Detection**
   - Make builds require manual configuration of MKL library paths and linking order
   - Windows MKL libraries are in different locations than Linux (`lib/` vs `lib/intel64/`)
   - CMake's `FindBLAS` and `FindLAPACK` modules automatically detect and configure MKL correctly

5. **MPI Integration**
   - Make builds require manual MPI compiler detection and flag configuration
   - CMake's `FindMPI` module automatically detects Intel MPI and configures it properly

### Advantages of CMake

- **Cross-platform**: Same build system works on Windows, Linux, and macOS
- **Automatic dependency detection**: Finds MKL, MPI, and other libraries automatically
- **Proper flag handling**: Separates compiler and linker flags correctly
- **Modern build system**: Better integration with IDE tools and CI/CD pipelines
- **Consistency**: Matches the build approach used by Quantum ESPRESSO
- **Maintainability**: Less platform-specific code to maintain

### Migration Path

The Make-based build (`make.inc`) was initially attempted but proved problematic on Windows with Intel oneAPI. The CMake build system, available in Wannier90's `develop` branch, provides a more robust and maintainable solution.

## Build Output

After a successful build, executables are located in:
```
upstream/wannier90/build-win-oneapi/
```
or
```
upstream/wannier90/build-win-oneapi/bin/
```

You should find:
- `wannier90.x` - Main Wannier90 executable
- `postw90.x` - Post-processing executable
- `w90chk2chk.x` - Checkpoint conversion utility (if tests are enabled)
- `w90spn2spn.x` - Spin conversion utility (if tests are enabled)

## Configuration

The CMake build is configured for:
- **Compiler**: Intel Fortran (ifx) and Intel C (icx)
- **Optimization**: Release mode (`/O2`)
- **Libraries**: Intel MKL (BLAS/LAPACK) with Intel threading
  - `mkl_intel_lp64` - 64-bit integer interface
  - `mkl_intel_thread` - Intel threading
  - `mkl_core` - Core MKL library
- **MPI**: Enabled by default (use `-NoMpi` flag for serial build)

The build script automatically:
- Sets `MKLROOT` and `CMAKE_PREFIX_PATH` from oneAPI installation
- Initializes oneAPI environment via `setvars.bat`
- Initializes Visual Studio environment (if found)
- Configures CMake with proper compiler paths

## Troubleshooting

### "Wannier90 source not found"

**Solution:** Run `.\scripts\refresh_wannier90_source.ps1` first to clone Wannier90.

### "setvars.bat not found"

**Cause:** Intel oneAPI not installed or not at default location.

**Solution:** 
1. Install Intel oneAPI HPC Toolkit
2. Or set `ONEAPI_ROOT` environment variable to your oneAPI installation path

### "MKL not found"

**Cause:** Intel MKL not installed or MKLROOT not set.

**Solution:**
1. Ensure oneAPI MKL component is installed (`intel.oneapi.win.mkl.devel`)
2. The build script should auto-detect MKL, but you can set `MKLROOT` manually:
   ```powershell
   $env:MKLROOT = "C:\Program Files (x86)\Intel\oneAPI\mkl\latest"
   ```

### "CMake not found" or "Ninja not found"

**Solution:** Install CMake and Ninja:
- **CMake**: Download from https://cmake.org/download/ or `choco install cmake`
- **Ninja**: `choco install ninja` or download from https://github.com/ninja-build/ninja/releases

### Build fails with linker errors

**Cause:** Visual Studio environment not properly initialized.

**Solution:** The script should automatically initialize VS. If it fails:
1. Install Visual Studio 2022 Build Tools with "Desktop development with C++" workload
2. Or manually run:
   ```cmd
   "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\Tools\VsDevCmd.bat" -arch=x64
   ```
   Then run the build script again.

## Notes

- **No patches required**: Unlike the MinGW build, this oneAPI build does not require any patches to Wannier90 source.
- **MPI build**: The default configuration builds with MPI support. Use `-NoMpi` flag to build a serial (non-MPI) version.
- **Executable naming**: Executables are built with the `.x` suffix (e.g., `wannier90.x`, `postw90.x`).
- **CMake build**: This build uses CMake, which is more modern, better integrated with Intel oneAPI, and follows the same approach as Quantum ESPRESSO. It requires the `develop` branch which includes CMake support.
- **Legacy Make build**: A Make-based build was attempted but encountered linker flag conflicts and Windows path issues. The CMake build system resolves these problems and is the recommended approach.
