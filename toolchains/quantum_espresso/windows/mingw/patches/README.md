# QE Windows MinGW Patches

This directory contains patches required to build Quantum ESPRESSO on Windows using MinGW-w64.

## Source

**All patches in this directory are sourced from the MSYS2 MINGW-packages repository:**

**Source Repository**: [https://github.com/msys2/MINGW-packages/tree/master/mingw-w64-quantum-espresso](https://github.com/msys2/MINGW-packages/tree/master/mingw-w64-quantum-espresso)

These patches are maintained by the MSYS2 project and are used in their official Quantum ESPRESSO package for MinGW-w64.

## Patches

### 001-use-srand-instead-of-srandom.patch
Replaces `srandom()` calls with `srand()` for MinGW compatibility. The `srandom()` function is not available on Windows/MinGW.

### 002-fix-build-devxlib-on-mingw.patch
Fixes compilation issues in devxlib for MinGW. **Note**: This patch must be applied **after** CMake configuration because `external/devxlib/src/timer.c` is downloaded during CMake config when `-DQE_DEVICEXLIB_INTERNAL=ON` is set.

The workflow applies a custom fix equivalent to this patch after CMake configuration.

### 003-fix-install-c-libraries.patch
Fixes installation of C libraries on Windows/MinGW.

### 004-fix-executable-suffix-on-mingw.patch
Fixes executable suffix handling in CMakeLists.txt files for Windows. Removes the `.x` suffix on Windows executables.

**Note**: This patch is applied in two stages:
1. Before CMake configuration (for files in the source tree)
2. After CMake configuration (for `external/pw2qmcpack/CMakeLists.txt` which is downloaded during CMake config)

### 005-fix-xmltools-eof-rewind.patch
Defensive fix for XML parsing: adds a rewind when EOF is reached to prevent "Read past ENDFILE record" crash. This prevents the hard crash but does not fix the root cause (empty tag being passed to `xmlr_opentag`).

### 006-fix-empty-tag-in-read-pp-header-v2.patch
**Root cause fix** for UPF v2.0.1 parsing crash on MinGW gfortran builds. Fixes two issues:

1. **`capitalize_if_v2` function**: Changed from using allocatable character concatenation in a loop (`strout = strout // ...`) to explicit allocation and direct character assignment. This avoids a compiler-specific bug in MinGW gfortran where allocatable character concatenation in loops can result in empty strings.

2. **`read_pp_header_v2` subroutine**: Stores the result of `capitalize_if_v2('pp_header')` in a local variable before passing it to `xmlr_readtag`, avoiding potential allocatable-to-fixed-length character conversion issues.

**Why MinGW triggered this but Intel oneAPI didn't:**
- **MinGW gfortran**: Has a known issue with allocatable character variables when using repeated concatenation in loops (`strout = strout // char`). The allocatable character may not be properly reallocated on each iteration, leading to an empty or corrupted result. Additionally, passing allocatable characters directly as arguments to procedures expecting fixed-length characters can cause conversion issues.
- **Intel oneAPI ifort**: More robust handling of allocatable character operations and better automatic conversion when passing allocatable characters to fixed-length character dummy arguments. The concatenation loop works correctly, and the direct function call in the argument list is handled properly.

The fix ensures compatibility with both compilers by using explicit allocation and direct character assignment, which is the most portable approach for allocatable character handling in Fortran.

## Usage

These patches are automatically applied by the CI workflow (`.github/workflows/qe-windows-mingw.yml`). They are applied in order using the `patch` command with `-Nbp1` flags:
- `-N`: Don't apply if already applied
- `-b`: Backup original files
- `-p1`: Strip first directory level from patch paths

## Updating Patches

To update patches from the MSYS2 repository:

```bash
cd toolchains/quantum_espresso/windows/mingw/patches
curl -L -o 001-use-srand-instead-of-srandom.patch "https://raw.githubusercontent.com/msys2/MINGW-packages/master/mingw-w64-quantum-espresso/001-use-srand-instead-of-srandom.patch"
curl -L -o 002-fix-build-devxlib-on-mingw.patch "https://raw.githubusercontent.com/msys2/MINGW-packages/master/mingw-w64-quantum-espresso/002-fix-build-devxlib-on-mingw.patch"
curl -L -o 003-fix-install-c-libraries.patch "https://raw.githubusercontent.com/msys2/MINGW-packages/master/mingw-w64-quantum-espresso/003-fix-install-c-libraries.patch"
curl -L -o 004-fix-executable-suffix-on-mingw.patch "https://raw.githubusercontent.com/msys2/MINGW-packages/master/mingw-w64-quantum-espresso/004-fix-executable-suffix-on-mingw.patch"
```

