# QE Windows MinGW Patches

This directory contains patches required to build Quantum ESPRESSO on Windows using MinGW-w64.

## For General Users

**Only two patches are necessary for building QE with MinGW:**

1. **002-fix-build-devxlib-on-mingw.patch** - Fixes compilation issues in devxlib
   - **Important**: This patch must be applied **after** CMake configuration and initial compilation, because `external/devxlib/src/timer.c` is downloaded during CMake config when `-DQE_DEVICEXLIB_INTERNAL=ON` is set. If you try to apply it before compilation, you may see "no destination file" errors. Run the build first, then apply this patch.

2. **005-mingw-upf-xml-fixes.patch** - Fixes UPF v2 XML parsing issues on MinGW gfortran
   - This patch can be applied before compilation.

**Note**: Patch 004 (executable suffix) is **not necessary** because the prefix change is already handled in the staging script.

## Source

**All patches in this directory are sourced from the MSYS2 MINGW-packages repository:**

**Source Repository**: [https://github.com/msys2/MINGW-packages/tree/master/mingw-w64-quantum-espresso](https://github.com/msys2/MINGW-packages/tree/master/mingw-w64-quantum-espresso)

These patches are maintained by the MSYS2 project and are used in their official Quantum ESPRESSO package for MinGW-w64.

## Patches

### Required Patches

#### 002-fix-build-devxlib-on-mingw.patch
Fixes compilation issues in devxlib for MinGW.

**⚠️ Important Timing**: This patch must be applied **after** CMake configuration and initial compilation attempt. The file `external/devxlib/src/timer.c` is downloaded during CMake config when `-DQE_DEVICEXLIB_INTERNAL=ON` is set. If you try to apply this patch before compilation, you will see "no destination file" errors because the file doesn't exist yet.

**Workflow:**
1. Run CMake configuration
2. Start compilation (this downloads devxlib)
3. Apply patch 002
4. Continue compilation

#### 005-mingw-upf-xml-fixes.patch
Comprehensive fix for UPF v2 XML parsing issues on MinGW gfortran builds. This unified patch consolidates all fixes from the deprecated patches 005-007.

**Fixes included:**
- Empty tag guard in `xmlr_opentag()` to prevent "match any '<'" bug
- `i2c()` function: changed from allocatable to fixed-length CHARACTER to avoid empty string bugs
- `capitalize_if_v2()` replaced with `capitalize_if_v2_into()` subroutine using fixed-length buffers
- All tag construction uses fixed-length buffers instead of allocatable characters
- Defensive EOF rewind to prevent "Read past ENDFILE record" crash
- Improved error handling in attribute parsing

**Why MinGW triggered this but Intel oneAPI didn't:**
- **MinGW gfortran**: Has a known issue with allocatable character variables when using repeated concatenation in loops (`strout = strout // char`). The allocatable character may not be properly reallocated on each iteration, leading to an empty or corrupted result. Additionally, passing allocatable characters directly as arguments to procedures expecting fixed-length characters can cause conversion issues.
- **Intel oneAPI ifort**: More robust handling of allocatable character operations and better automatic conversion when passing allocatable characters to fixed-length character dummy arguments.

### Optional/Deprecated Patches

#### 001-use-srand-instead-of-srandom.patch
Replaces `srandom()` calls with `srand()` for MinGW compatibility. The `srandom()` function is not available on Windows/MinGW.

**Status**: May be needed depending on QE version, but not currently required for QE 7.5.

#### 003-fix-install-c-libraries.patch
Fixes installation of C libraries on Windows/MinGW.

**Status**: May be needed depending on build configuration.

#### 004-fix-executable-suffix-on-mingw.patch
Fixes executable suffix handling in CMakeLists.txt files for Windows. Removes the `.x` suffix on Windows executables.

**Status**: **NOT NECESSARY** - The prefix change is already handled in the staging script, so this patch is redundant.

## Usage

### Manual Application

Apply patches using the `patch` command with `-Nbp1` flags:
- `-N`: Don't apply if already applied
- `-b`: Backup original files
- `-p1`: Strip first directory level from patch paths

**Example:**
```bash
# Apply patch 005 (can be done before compilation)
patch -Nbp1 < patches/005-mingw-upf-xml-fixes.patch

# After CMake config and initial compilation (when devxlib is downloaded)
patch -Nbp1 < patches/002-fix-build-devxlib-on-mingw.patch
```

### Automated Application

These patches are automatically applied by the CI workflow (`.github/workflows/qe-windows-mingw.yml`). The workflow handles the timing correctly:
- Patch 005 is applied before compilation
- Patch 002 is applied after CMake configuration (when devxlib is available)

## Updating Patches

To update patches from the MSYS2 repository:

```bash
cd toolchains/quantum_espresso/windows/mingw/patches
curl -L -o 002-fix-build-devxlib-on-mingw.patch "https://raw.githubusercontent.com/msys2/MINGW-packages/master/mingw-w64-quantum-espresso/002-fix-build-devxlib-on-mingw.patch"
```

**Note**: Patch 005 (`005-mingw-upf-xml-fixes.patch`) is a custom unified patch that consolidates fixes from multiple MSYS2 patches. It is maintained in this repository and should not be updated from MSYS2.

For optional patches (if needed for other QE versions):
```bash
curl -L -o 001-use-srand-instead-of-srandom.patch "https://raw.githubusercontent.com/msys2/MINGW-packages/master/mingw-w64-quantum-espresso/001-use-srand-instead-of-srandom.patch"
curl -L -o 003-fix-install-c-libraries.patch "https://raw.githubusercontent.com/msys2/MINGW-packages/master/mingw-w64-quantum-espresso/003-fix-install-c-libraries.patch"
```

