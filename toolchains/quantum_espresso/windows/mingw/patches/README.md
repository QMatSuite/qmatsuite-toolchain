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

