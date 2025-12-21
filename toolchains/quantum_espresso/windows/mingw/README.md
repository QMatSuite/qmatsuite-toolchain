# Quantum ESPRESSO on Windows (MSYS2 UCRT64)

Windows-specific build guide for Quantum ESPRESSO 7.5 using MSYS2 + MinGW-w64 (UCRT64 toolchain).

## Quick Start

**Goal:** Build and run QE in <10 minutes using known-good production flags.

### Prerequisites

1. **Install MSYS2** (if not already installed):
   - Download from https://www.msys2.org/
   - Or use Chocolatey: `choco install msys2`

2. **Open MSYS2 UCRT64 terminal** (not the regular MSYS2 terminal)

3. **Install required packages:**
   ```bash
   pacman -Syu
   pacman -S --needed \
     base-devel git wget patch make pkgconf \
     mingw-w64-ucrt-x86_64-toolchain \
     mingw-w64-ucrt-x86_64-gcc-fortran \
     mingw-w64-ucrt-x86_64-openblas \
     mingw-w64-ucrt-x86_64-fftw \
     mingw-w64-ucrt-x86_64-hdf5 \
     mingw-w64-ucrt-x86_64-pkgconf \
     mingw-w64-ucrt-x86_64-ntldd \
     mingw-w64-ucrt-x86_64-binutils
   ```

### Build Steps

```bash
# 1. Navigate to QE source directory
cd /c/Users/$(whoami)/qmatsuite-toolchain/toolchains/quantum_espresso/windows/mingw/upstream/q-e-qe-7.5

# 2. Apply required patches
PATCH_DIR="/c/Users/$(whoami)/qmatsuite-toolchain/toolchains/quantum_espresso/windows/mingw/patches"
patch -Nbp1 -i "$PATCH_DIR/001-use-srand-instead-of-srandom.patch"
patch -Nbp1 -i "$PATCH_DIR/002-fix-build-devxlib-on-mingw.patch"
patch -Nbp1 -i "$PATCH_DIR/004-fix-executable-suffix-on-mingw.patch"
patch -Nbp1 -i "$PATCH_DIR/005-fix-xmltools-eof-rewind.patch"
patch -Nbp1 -i "$PATCH_DIR/006-fix-empty-tag-robust-mingw.patch"
patch -Nbp1 -i "$PATCH_DIR/007-fix-i2c-mingw-empty-string.patch"

# 3. Set production compiler flags (PRODUCTION-SAFE)
export FFLAGS="-O2 -fopenmp -Wno-missing-include-dirs -w -std=legacy -fallow-argument-mismatch -fallow-invalid-boz -fmax-stack-var-size=0"
export CFLAGS="-O2 -Wno-incompatible-pointer-types"
export LDFLAGS="-Wl,--stack,268435456 -fopenmp"
export QE_PREFIX="/c/Users/$(whoami)/qmatsuite-toolchain/toolchains/quantum_espresso/windows/mingw/upstream/q-e-qe-7.5/artifacts"

# 4. Configure
./configure \
  MPIF90=gfortran \
  BLAS_LIBS="-lopenblas" \
  LAPACK_LIBS="-lopenblas" \
  FFT_LIBS="-lfftw3 -lfftw3_omp" \
  DFLAGS="-D_WIN32 -D__FFTW3" \
  --enable-openmp \
  --prefix="$QE_PREFIX"

# 5. Build (single-threaded to avoid race conditions)
make veryclean
make -j1

# 6. Install
make install

# 7. Stage executables and DLLs (from PowerShell)
cd /c/Users/$(whoami)/qmatsuite-toolchain/toolchains/quantum_espresso/windows/mingw
powershell -Command ".\scripts\stage_qe_dlls.ps1"

# 8. Test (critical: test with multiple threads!)
export OMP_NUM_THREADS=4
cd dist/win-mingw
./pw.exe -in scf-cg.in > scf-cg.out 2>&1
tail -20 scf-cg.out
```

**Important Notes:**
- **Always use UCRT64 packages** (`mingw-w64-ucrt-x86_64-*`), not mingw64
- **Stack reserve is critical**: `-Wl,--stack,268435456` (256MB) is required
- **Test with OMP_NUM_THREADS>1** before declaring success
- **Patches 005-007 are required** for MinGW robustness (UPF v2 parsing)

---

## Production Build Guide

### Prerequisites

#### 1. Install MSYS2

Install MSYS2 using Chocolatey:
```powershell
choco install msys2
```

Or download and install manually from: https://www.msys2.org/

#### 2. Install Required Packages

Open **MSYS2 UCRT64 terminal** (not the regular MSYS2 terminal) and run:

```bash
# Update package database
pacman -Syu

# Install base development tools (MSYS packages, no prefix)
pacman -S --needed base-devel git wget patch make pkgconf

# Install UCRT64 toolchain and libraries (all ucrt64-prefixed)
pacman -S --needed \
  mingw-w64-ucrt-x86_64-toolchain \
  mingw-w64-ucrt-x86_64-gcc-fortran \
  mingw-w64-ucrt-x86_64-openblas \
  mingw-w64-ucrt-x86_64-fftw \
  mingw-w64-ucrt-x86_64-hdf5 \
  mingw-w64-ucrt-x86_64-pkgconf \
  mingw-w64-ucrt-x86_64-ntldd \
  mingw-w64-ucrt-x86_64-binutils
```

**Critical:** Use **UCRT64** packages (`mingw-w64-ucrt-x86_64-*`), not mingw64. MSYS2 environments are separate and do not share packages.

### Step 1: Download and Extract QE Source

In the MSYS2 UCRT64 terminal:

```bash
# Create upstream directory (adjust path as needed)
cd /c/Users/$(whoami)/qmatsuite-toolchain/toolchains/quantum_espresso/windows/mingw
mkdir -p upstream
cd upstream

# Download QE tarball (if not already downloaded)
QE_VERSION="qe-7.5"
QE_URL="https://gitlab.com/QEF/q-e/-/archive/${QE_VERSION}/q-e-${QE_VERSION}.tar.gz"

if [ ! -f "qe.tar.gz" ]; then
    wget -q "$QE_URL" -O qe.tar.gz
fi

# Extract (if not already extracted)
if [ ! -d "q-e-qe-7.5" ]; then
    tar -xzf qe.tar.gz
    EXTRACTED=$(find . -maxdepth 1 -type d -name "q-e-${QE_VERSION}*" | head -1)
    if [ -n "$EXTRACTED" ]; then
        mv "$EXTRACTED" q-e-qe-7.5
    fi
fi

cd q-e-qe-7.5
```

### Step 2: Apply Required Patches

Apply the MinGW-specific patches from the repository:

```bash
# Navigate to QE source directory
cd /c/Users/$(whoami)/qmatsuite-toolchain/toolchains/quantum_espresso/windows/mingw/upstream/q-e-qe-7.5

# Set patch directory
PATCH_DIR="/c/Users/$(whoami)/qmatsuite-toolchain/toolchains/quantum_espresso/windows/mingw/patches"

# Apply patches (use -N to skip if already applied)
patch -Nbp1 -i "$PATCH_DIR/001-use-srand-instead-of-srandom.patch"
patch -Nbp1 -i "$PATCH_DIR/002-fix-build-devxlib-on-mingw.patch"
patch -Nbp1 -i "$PATCH_DIR/004-fix-executable-suffix-on-mingw.patch"

# Critical MinGW robustness fixes for UPF v2 parsing
patch -Nbp1 -i "$PATCH_DIR/005-fix-xmltools-eof-rewind.patch"
patch -Nbp1 -i "$PATCH_DIR/006-fix-empty-tag-robust-mingw.patch"
patch -Nbp1 -i "$PATCH_DIR/007-fix-i2c-mingw-empty-string.patch"
```

**Why patches 005-007 are required:** These fix MinGW gfortran-specific issues with allocatable CHARACTER functions that cause UPF v2 parsing to fail. See Appendix A for details.

### Step 3: Configure QE

Configure QE using the traditional `./configure` script with production-safe flags:

```bash
# Production-safe compiler flags (see "Compiler Flags Reference" below)
# Might be useful if problem with stack overflow: FFLAGS = -fmax-stack-var-size=0
export FFLAGS="-O2 -fopenmp -Wno-missing-include-dirs -w -std=legacy -fallow-argument-mismatch -fallow-invalid-boz "
export CFLAGS="-O2 -Wno-incompatible-pointer-types"
export LDFLAGS="-Wl,--stack,268435456 -fopenmp"
export QE_PREFIX="/c/Users/$(whoami)/qmatsuite-toolchain/toolchains/quantum_espresso/windows/mingw/upstream/q-e-qe-7.5/artifacts"

./configure \
  MPIF90=gfortran \
  BLAS_LIBS="-lopenblas" \
  LAPACK_LIBS="-lopenblas" \
  FFT_LIBS="-lfftw3 -lfftw3_omp" \
  DFLAGS="-D_WIN32 -D__FFTW3" \
  --enable-openmp \
  --prefix="$QE_PREFIX"
```

**Configuration options explained:**
- `MPIF90=gfortran`: Use gfortran for MPI (serial build)
- `BLAS_LIBS="-lopenblas"`: Link against OpenBLAS
- `LAPACK_LIBS="-lopenblas"`: Use OpenBLAS for LAPACK
- `FFT_LIBS="-lfftw3 -lfftw3_omp"`: Use external FFTW3 with OpenMP support
- `DFLAGS="-D_WIN32 -D__FFTW3"`: Windows-specific flags and external FFTW3
- `--enable-openmp`: Enable OpenMP support

**Important:** When passing `FFT_LIBS` to configure, you must also set `-D__FFTW3` in `DFLAGS` manually. Otherwise configure may fail to set the FFT backend macro, resulting in compilation errors.

### Step 4: Build QE

Build with single-threaded make to avoid race conditions:

```bash
# Clean previous build
make veryclean

# Build (single-threaded to avoid module file race conditions)
make -j1

# Or build specific target (e.g., pw.x)
make pw -j1
```

**Why `-j1`:** Fortran module files can have race conditions during parallel builds. Single-threaded builds are more reliable.

### Step 5: Install

Install to the prefix directory:

```bash
make install
```

Executables will be installed to `$QE_PREFIX/bin/` (e.g., `artifacts/bin/pw.x`).

### Step 6: Stage Executables and Required DLLs

After building and installing, stage the executables and required DLLs to a clean distribution directory using the staging script:

**From PowerShell:**
```powershell
cd toolchains\quantum_espresso\windows\mingw
.\scripts\stage_qe_dlls.ps1
```

**Or with custom paths:**
```powershell
.\scripts\stage_qe_dlls.ps1 -ArtifactsDir "C:\path\to\qe\artifacts" -DistDir "C:\path\to\dist" -MingwBinDir "C:\msys64\ucrt64\bin"
```

**What the staging script does:**

1. **Finds executables**: Recursively searches for all `*.x` files in the artifacts directory
2. **Copies to distribution directory**: Copies all `*.x` files to `dist\win-mingw` in a flat layout
3. **Renames to .exe**: Renames all copied `*.x` files to `*.exe` (e.g., `pw.x` → `pw.exe`)
4. **Copies required DLLs**: Copies required runtime DLLs from MSYS2 UCRT64 bin directory

**MSYS2 UCRT64 directory detection:**

The script automatically finds your MSYS2 UCRT64 installation in this order:
- `C:\msys64\ucrt64\bin`
- `C:\tools\msys64\ucrt64\bin`
- `$env:MSYS2_ROOT\ucrt64\bin` (if `MSYS2_ROOT` is set)
- `%ProgramFiles%\msys64\ucrt64\bin`

**Important:** The script uses **UCRT64** bin directory, not mingw64. This matches the UCRT64 toolchain used for building.

**Staged DLLs:**

The script copies these DLL patterns from UCRT64 bin:
- `libgcc_s_seh-1.dll` - GCC runtime
- `libwinpthread-1.dll` - POSIX threads
- `libgfortran-*.dll` - GNU Fortran runtime
- `libquadmath-*.dll` - Quadmath library
- `libgomp-1.dll` - OpenMP runtime
- `libstdc++-6.dll` - C++ standard library
- `libopenblas*.dll` - OpenBLAS runtime
- `libfftw3-3.dll` - FFTW3 runtime
- `libfftw3_omp-3.dll` - FFTW3 OpenMP support

**Distribution directory:**

The staging script creates a clean distribution directory at:
```
toolchains\quantum_espresso\windows\mingw\dist\win-mingw
```

This directory contains only `.exe` and `.dll` files - all executables and their required runtime DLLs in a flat layout, ready for distribution or testing.

**Note on ntldd:**

The staging script uses a fixed, practical DLL list rather than automatic detection via `ntldd`. If you need `ntldd` for other purposes, install it with:
```bash
pacman -S mingw-w64-ucrt-x86_64-ntldd
```

### Step 7: Verify Build

**Verify stack reserve size:**
```bash
objdump -x artifacts/bin/pw.x | grep -i SizeOfStackReserve
# Should show: SizeOfStackReserve 0000000010000000 (256MB)
```

**Verify which binary you're running:**
```bash
which pw.exe
ls -lh "$(which pw.exe)"
```

**Test with multiple threads (critical!):**
```bash
export OMP_NUM_THREADS=4
cd dist/win-mingw
./pw.exe -in scf-cg.in > scf-cg.out 2>&1
tail -20 scf-cg.out
```

---

## Compiler Flags Reference

### Production-Safe Flags

These flags are safe for production builds and OpenMP execution:

**FFLAGS (Production):**
```bash
-O2 -fopenmp -Wno-missing-include-dirs -w -std=legacy -fallow-argument-mismatch -fallow-invalid-boz -fmax-stack-var-size=0
```

- `-O2`: Optimization level 2
- `-fopenmp`: Enable OpenMP support
- `-Wno-missing-include-dirs`: Suppress missing include directory warnings
- `-w`: Suppress all warnings
- `-std=legacy`: Use legacy Fortran standard (for compatibility)
- `-fallow-argument-mismatch`: Allow BLAS argument type mismatches
- `-fallow-invalid-boz`: Allow invalid binary/octal/hexadecimal constants
- `-fmax-stack-var-size=0`: Push large automatic arrays off stack (Windows-specific)

**CFLAGS (Production):**
```bash
-O2 -Wno-incompatible-pointer-types
```

- `-O2`: Optimization level 2
- `-Wno-incompatible-pointer-types`: Suppress incompatible pointer type warnings

**LDFLAGS (Production):**
```bash
-Wl,--stack,268435456 -fopenmp
```

- `-Wl,--stack,268435456`: Set stack reserve to 256MB (critical for Windows)
- `-fopenmp`: Link OpenMP runtime

### Debug-Only Flags (NOT OpenMP Safe)

**⚠️ WARNING:** These flags are useful for single-threaded debugging but can destabilize OpenMP execution with multiple threads. Do NOT use in production builds.

**Debug flags (single-thread only):**
```bash
-g -fbacktrace -fcheck=all -finit-character=0 -finit-local-zero -ffpe-trap=invalid,zero,overflow
```

- `-g`: Include debug symbols
- `-fbacktrace`: Enable backtrace on runtime errors
- `-fcheck=all`: Enable all runtime checks (bounds, pointers, memory, etc.)
- `-finit-character=0`: Initialize CHARACTER variables to null
- `-finit-local-zero`: Initialize local variables to zero
- `-ffpe-trap=invalid,zero,overflow`: Trap floating-point exceptions

**Known issue with debug flags:**

Using `-fcheck=all` or `-ffpe-trap` with `OMP_NUM_THREADS>1` can cause false "recursive call" errors in XC routines (e.g., `slater`, `pz`). These routines are not recursive; the error is a false positive from runtime checks under OpenMP.

**Recommendation:** Use debug flags only for single-threaded diagnosis. Always validate production builds with `OMP_NUM_THREADS>1`.

---

## Running and Testing QE

### From PowerShell

```powershell
$env:PATH="$HOME\qmatsuite-toolchain\toolchains\quantum_espresso\windows\mingw\dist\win-mingw;$env:PATH"
pw.exe -in si100.in > "si100_mingw_omp1.out" 2>&1
```

### From MSYS2 UCRT64 Terminal

```bash
# Navigate to distribution directory
cd /c/Users/$(whoami)/qmatsuite-toolchain/toolchains/quantum_espresso/windows/mingw/dist/win-mingw

# Run QE
./pw.exe -in input.in > output.out 2>&1
```

### From Windows Command Prompt

1. Navigate to the distribution directory:
   ```cmd
   cd C:\Users\%USERNAME%\qmatsuite-toolchain\toolchains\quantum_espresso\windows\mingw\dist\win-mingw
   ```

2. Run QE:
   ```cmd
   pw.exe -in input.in > output.out 2>&1
   ```

**Note:** All required DLLs are in the same directory as the executables, so no PATH modification is needed.

### Testing with Multiple Threads

**Critical:** Always test with `OMP_NUM_THREADS>1` to verify OpenMP stability:

```bash
export OMP_NUM_THREADS=4
pw.exe -in scf-cg.in > scf-cg.out 2>&1
```

If the build works with `OMP_NUM_THREADS=1` but fails with `OMP_NUM_THREADS>=2`, check for debug flags in your build configuration.

---

## Required DLLs

The following DLLs are automatically staged by `stage_qe_dlls.ps1` from MSYS2 UCRT64:

**Runtime DLLs (from UCRT64 bin):**
- `libgcc_s_seh-1.dll` - GCC runtime library
- `libwinpthread-1.dll` - POSIX threads implementation
- `libgfortran-*.dll` - GNU Fortran runtime (typically `libgfortran-5.dll`)
- `libquadmath-*.dll` - Quadmath library (required by libgfortran)
- `libgomp-1.dll` - OpenMP runtime (required if OpenMP enabled)
- `libstdc++-6.dll` - C++ standard library
- `libopenblas*.dll` - OpenBLAS runtime (BLAS/LAPACK implementation)
- `libfftw3-3.dll` - FFTW3 runtime
- `libfftw3_omp-3.dll` - FFTW3 OpenMP support

**Important:** These DLLs must be from the **UCRT64** toolchain (`C:\msys64\ucrt64\bin`), not mingw64, to match the build toolchain.

**Note:** The `stage_qe_dlls.ps1` script automatically copies all required DLLs to the distribution directory. No manual copying is needed.

---

## Troubleshooting

### Build Errors

1. **Compiler not found**: Ensure you're using MSYS2 UCRT64 terminal, not the regular MSYS2 terminal.

2. **Library not found**: Check that all required packages are installed:
   ```bash
   pacman -Qs mingw-w64-ucrt-x86_64-openblas
   pacman -Qs mingw-w64-ucrt-x86_64-fftw
   ```

3. **Patch application fails**: Patches may already be applied. Use `-N` flag to skip if already applied.

### Common Compilation Failures

#### 1. FFTXlib / fft_scalar.f90: No FFT backend selected

**Error:**
```
fft_scalar.f90: fatal error: #error No fft_scalar backend selected!
```

**Cause:** No FFT backend macro defined (`-D__FFTW3` / `-DFFTW3` / `-DFFTW` / `-DDFTI`) because `configure` didn't detect FFTW3 libraries/includes.

**Fix:**
```bash
# Ensure FFTW3 is installed in UCRT64
pacman -S --needed mingw-w64-ucrt-x86_64-fftw

# Verify FFTW3 is available
pkg-config --libs fftw3

# Reconfigure with explicit FFT flags and DFLAGS
export FFT_LIBS="-lfftw3 -lfftw3_omp"
export DFLAGS="-D_WIN32 -D__FFTW3"
./configure ... (other options) ...
```

**Verification:**
```bash
# Check that make.inc contains exactly one FFT backend macro
grep -E "DFLAGS.*FFTW|DFLAGS.*DFTI" make.inc
# Should show: -D__FFTW3 (or -DFFTW3, -DFFTW, -DDFTI) - exactly one!
```

#### 2. LAXlib / ptoolkit.f90: Type mismatch in zgemm call

**Error:**
```
ptoolkit.f90: error: Type mismatch in argument 'c' at (1) (COMPLEX(8)/REAL(8))
```

**Cause:** Newer gfortran treats legacy BLAS argument mismatches as hard errors.

**Fix:**
```bash
# Add -fallow-argument-mismatch to Fortran flags
export FFLAGS="-O2 -fopenmp -Wno-missing-include-dirs -w -std=legacy -fallow-argument-mismatch -fallow-invalid-boz -fmax-stack-var-size=0"

# Reconfigure and rebuild
make veryclean
./configure ... (your options) ...
make -j1
```

**Verification:**
```bash
# Check that FFLAGS contains the flag
grep "fallow-argument-mismatch" make.inc
```

#### 3. PW/src: Cannot open module file 'device_fbuff_m.mod'

**Error:**
```
Fatal Error: Cannot open module file 'device_fbuff_m.mod' for reading
```

**Cause:** Fortran module `.mod` file not generated/visible due to missing `-J` module output directory flag and/or parallel build race condition.

**Fix:**
```bash
# Add -J flag to specify module output directory
export F90FLAGS="-J$QE_PREFIX/include"
export FFLAGS="$FFLAGS -I$QE_PREFIX/include"

# Clean rebuild (use -j1 to avoid race conditions)
make veryclean
./configure ... (your options) ...
make -j1  # Single-threaded build to avoid race conditions
```

**Verification:**
```bash
# Check that module file exists
find . -name "device_fbuff_m.mod"
# Should show: ./include/device_fbuff_m.mod (or similar)

# Check that make.inc contains -J flag
grep -E "F90FLAGS.*-J|FFLAGS.*-J" make.inc
```

**Workaround:** If the issue persists, build with `make -j1` (single-threaded) to confirm it's a race condition, then investigate the build order.

### Runtime Errors

1. **DLL not found**: Ensure all required DLLs are in the same directory as the executables. The staging script should have copied them automatically.

2. **Path issues**: Use forward slashes in MSYS2 terminal, or use Windows paths with proper escaping.

3. **Stack overflow (0xC00000FD)**: Verify stack reserve size:
   ```bash
   objdump -x pw.exe | grep -i SizeOfStackReserve
   # Should show: SizeOfStackReserve 0000000010000000 (256MB)
   ```
   If it shows `0000000000200000` (2MB), rebuild with `-Wl,--stack,268435456` in LDFLAGS.

4. **Wrong binary being run**: Always verify which binary you're executing:
   ```bash
   which pw.exe
   ls -lh "$(which pw.exe)"
   ```

---

## Known Issues

### UPF/XML Parsing Crash - RESOLVED

**Status:** Fixed via patches 005, 006, and 007. See Appendix A for full details.

The issue was caused by MinGW gfortran mishandling allocatable CHARACTER function returns in tag construction. The fix replaces allocatable functions with fixed-length buffer subroutines.

**Required patches:**
- `005-fix-xmltools-eof-rewind.patch` — defensive EOF rewind fix
- `006-fix-empty-tag-robust-mingw.patch` — `capitalize_if_v2` fix
- `007-fix-i2c-mingw-empty-string.patch` — `i2c()` fix and tag construction updates

**Symptom (before fix):**

On Windows MSYS2 UCRT64 (gfortran), both the MSYS2 pacman-installed `pw.exe` and locally compiled `pw.x` would fail when reading UPF v2 pseudopotentials with:

```
end of file reached, tag  not found
end of file reached, tag not found
At line 1081 of file xmltools.f90
Fortran runtime error: Read past ENDFILE record
```

**Solution:** Apply patches 005-007 before building. These patches are now part of the standard build process.

---

## Patches

The build process uses patches to fix MinGW-specific build issues. These patches are stored locally in `toolchains/quantum_espresso/windows/mingw/patches/`.

**Required patches:**
- `001-use-srand-instead-of-srandom.patch` - Use `srand` instead of `srandom` (Windows compatibility)
- `002-fix-build-devxlib-on-mingw.patch` - Fix devxlib build on MinGW
- `004-fix-executable-suffix-on-mingw.patch` - Fix executable suffix handling
- `005-fix-xmltools-eof-rewind.patch` - Defensive EOF rewind fix (UPF parsing)
- `006-fix-empty-tag-robust-mingw.patch` - Fix allocatable CHARACTER bug (UPF parsing)
- `007-fix-i2c-mingw-empty-string.patch` - Fix `i2c()` allocatable CHARACTER bug (UPF parsing)

See `toolchains/quantum_espresso/windows/mingw/patches/README.md` for details about each patch.

---

## CI Build (CMake)

The CI workflow (`.github/workflows/qe-windows-mingw.yml`) uses CMake instead of the configure script. For CMake-based builds, see the workflow file for the exact configuration options used.

---

## Migration from Intel oneAPI

This toolchain replaces the previous Intel oneAPI/ifort/MS-MPI based Windows builds. The MinGW-w64 toolchain provides:

- Open-source compiler (gfortran)
- No proprietary dependencies
- Better compatibility with open-source scientific software
- Easier redistribution (no Intel license restrictions)

---

## Appendix A: Complete Failure Timeline and Postmortem

This appendix documents the complete journey from initial build attempts to a working production configuration, including all major failures, diagnostic steps, and fixes. This is essential reading for understanding why each step exists and for debugging similar issues.

### Phase 0: Initial Symptom - Runtime "Silent Crash"

**Symptom:**

Both locally built QE and the official MSYS2-packaged QE binaries (`mingw-w64-ucrt-x86_64-quantum-espresso`) would appear to "hang" or crash immediately after startup with no useful error output. The process would print the QE banner and "Current dimensions of program PWSCF are: ..." then stop with no further output.

**Why this matters:**

Because this also happened with the MSYS2 official build, we concluded it was not just our local build mistakes, and suspected a Windows stack/runtime configuration issue (similar to a previous oneAPI case).

**Evidence:**

- Runtime directly dies / nothing printed for some `pw.exe` runs
- Process terminates with Windows exception `0xC00000FD` (STATUS_STACK_OVERFLOW)
- Verified we were sometimes accidentally running an old `pw.exe` from `/ucrt64/bin` rather than our freshly built `pw.x`

**Must-do sanity checks:**

```bash
# Always verify which binary you're running
which pw.exe
ls -lh "$(which pw.exe)"

# Check stack reserve size (requires binutils)
pacman -S --needed mingw-w64-ucrt-x86_64-binutils
objdump -x "$(which pw.exe)" | grep -i SizeOfStackReserve
objdump -x pw.x | grep -i SizeOfStackReserve
```

**Key observation:**

```
/ucrt64/bin/pw.exe
SizeOfStackReserve 0000000000200000   (2MB - too small!)

pw.x (our build)
SizeOfStackReserve 0000000010000000   (256MB - correct)
```

**Conclusion:**

On Windows, stack reserve can differ between binaries and can be a root cause of "instant death". The default 2MB stack reserve is insufficient for QE's large automatic arrays. Make sure we always run the intended binary, not the stale system one.

---

### Phase 1: Stack-Related Fixes to Get Past Early Runtime Death

**Context:**

We had a previous oneAPI build that required stack/noheap adjustments. For gfortran we don't have oneAPI's "noheap", so we used the link-time stack reserve approach.

**Solution:**

Set link flags to increase stack reserve:

```bash
export LDFLAGS="-Wl,--stack,268435456 -fopenmp"
```

Also note the OpenMP stack variable behavior:

```bash
export FFLAGS="-O2 -fopenmp -fmax-stack-var-size=0 ..."
export FCFLAGS="-O2 -fopenmp -fmax-stack-var-size=0"
```

**Why this matters:**

On Windows, large automatic arrays can trigger stack overflow or undefined behavior. `-Wl,--stack,<bytes>` increases reserve stack size (256MB = 268435456 bytes). `-fmax-stack-var-size=0` pushes big locals off stack (compiler-dependent behavior, but helps).

**Result:**

This phase got us to the point where QE started printing banners / progressed further instead of dying instantly.

---

### Phase 2: FFTW/OpenMP Link Failures

**Symptom:**

After getting the build further, we hit link errors involving FFTW threads.

**Error message (verbatim):**

```
undefined reference to `fftw_init_threads'
undefined reference to `fftw_plan_with_nthreads'
collect2.exe: error: ld returned 1 exit status
```

**Link line context:**

The configure script tries to link with:
```
... -lfftw3 -lfftw3_omp ...
```

**Root cause:**

FFTW3 was not installed in UCRT64. It had been installed previously in another MSYS2 environment (mingw64), and we assumed it carried over.

**Key lesson:**

MSYS2 environments do not share packages (mingw64 vs ucrt64 are separate). Installing FFTW in mingw64 does NOT install it in ucrt64.

**How to verify FFTW presence:**

```bash
# Check if FFTW3 is available
pkg-config --libs fftw3
pkg-config --variable=prefix fftw3

# After building, verify DLL dependencies
ntldd -R pw.x | head -n 40
```

**Solution:**

Install FFTW in UCRT64 (use ucrt64-prefixed packages consistently):

```bash
pacman -S --needed mingw-w64-ucrt-x86_64-fftw
```

**Result:**

Build completes successfully after installing FFTW3 in the correct MSYS2 environment.

---

### Phase 3: Configure Gotcha - --enable-openmp Typo

**Symptom:**

Configure step fails with a confusing error.

**Error message (verbatim):**

```
-bash: -enable-openmp: command not found
```

**Root cause:**

We mistakenly wrote `-enable-openmp` (single dash) instead of `--enable-openmp` (double dash). Bash treated it as a command.

**Solution:**

Use double dash for configure options:

```bash
./configure --enable-openmp ...
```

**Why this matters:**

QE configure option is `--enable-openmp` (double dash). This is a common mistake when copying configure commands.

---

### Phase 4: The Real Major Blocker - UPF v2 XML Parsing Crash on MinGW

**Symptom:**

Even after FFTW/OpenMP issues were fixed, runtime failed during reading UPF v2 files.

**Original runtime failure signature (verbatim):**

```
end of file reached, tag  not found
end of file reached, tag not found
At line 1081 of file xmltools.f90
Fortran runtime error: Read past ENDFILE record
```

**Initial investigation:**

We confirmed the UPF file existed and worked with oneAPI / other binaries, so the file was not the problem:
- It would not say "file not found"
- Absolute path didn't help
- `dos2unix` didn't help

**Code inspection:**

File: `upflib/xmltools.f90`

Line 1081:
```fortran
read(xmlunit,'(a)', end=10) line
```

**Key reasoning:**

The error message had **two spaces**: `tag  not found`, implying `trim(tag)` was empty or whitespace. This indicates:
- The XML reader was being asked to open/read an **empty tag name**
- Not that the UPF XML content was missing

**Diagnostic guards added:**

We added diagnostics/guards and captured a backtrace proving the call chain:

```
FATAL: xmlr_opentag called with empty/whitespace tag
Tag value (repr): []
xmlunit=-10, nlevel=1
open tags: [root] [UPF]
```

Backtrace showed:
```
xmlr_opentag -> readtag_c -> read_pp_header_v2 -> read_upf_new -> readpp -> pwscf
```

**Root cause hypothesis:**

MinGW gfortran mishandles some combinations of:
- allocatable CHARACTER function returns
- loop concatenation (`strout = strout // char`)
- passing function results directly as actual arguments

This produced empty strings for tag names like `'pp_header'` in UPF v2 mode.

**The fix implemented (canonical solution):**

We replaced the allocatable-return function with a fixed-length buffer subroutine and updated all call sites.

**Conceptually:**

Replace:
```fortran
tag = capitalize_if_v2('pp_header', v2)
call xmlr_readtag(tag, ...)
```

With:
```fortran
CHARACTER(LEN=64) :: tag_buf
CALL capitalize_if_v2_into('pp_header', v2, tag_buf)
CALL xmlr_readtag(TRIM(tag_buf), ...)
```

The "into" subroutine uses direct character assignment (no concatenation):
```fortran
DO n = 1, lt
  strout(n:n) = capital(strin(n:n))
END DO
```

**Files touched:**

- `upflib/read_upf_new.f90` — updated all call sites (27 instances)
- `upflib/write_upf_new.f90` — updated all call sites for consistency (26 instances)
- `upflib/read_upf_new.f90` (the capitalization helper implementation)
- `upflib/upf_utils.f90` provides `capital()` (used by helper)

**Additional fix for i2c() function:**

The `i2c()` function in `xmltools.f90` also returned allocatable CHARACTER, which caused empty strings when constructing tags like `'PP_BETA.'//i2c(nb)`. We fixed this by:
- Changing `i2c()` to return fixed-length `CHARACTER(LEN=32)` instead of allocatable
- Replacing all `i2c()` usages in tag construction with direct `WRITE` to fixed buffers

**Patches applied:**

- `005-fix-xmltools-eof-rewind.patch` — defensive EOF rewind fix
- `006-fix-empty-tag-robust-mingw.patch` — `capitalize_if_v2` fix
- `007-fix-i2c-mingw-empty-string.patch` — `i2c()` fix and tag construction updates

**Why this matters:**

This is a MinGW-specific robustness fix; oneAPI didn't fail here. We keep it because it makes the code more portable and avoids a real compiler edge case.

**Result:**

UPF v2 files now parse correctly. QE progresses past pseudopotential reading.

---

### Phase 5: After XML Fix - New Runtime Issue Only When OMP Threads >= 2

**Symptom:**

Once UPF parsing was fixed, runtime would proceed much further and read pseudopotentials correctly, but crash when OpenMP used more than 1 thread.

**Evidence:**

- `OMP_NUM_THREADS=1` → test passes
- `OMP_NUM_THREADS>=2` → runtime error around XC routines

**Representative error (verbatim):**

```
Fortran runtime error: Recursive call to nonrecursive procedure 'slater'
At line 16 of file qe_funct_exch_lda_lsda.f90
```

And similarly for `pz` in correlation.

**Source code verification:**

File: `XClib/qe_funct_exch_lda_lsda.f90`

```fortran
SUBROUTINE slater( rs, ex, vx )
  ...
  ex = f * alpha / rs
  vx = 4._DP / 3._DP * f * alpha / rs
END SUBROUTINE slater
```

The procedure is **not recursive**.

**Why this happened:**

This was caused by debugging/runtime-check flags we had added earlier while diagnosing the XML bug, which were not OpenMP-safe in this context (or triggered gfortran runtime behavior under multi-threading).

**Problematic debug-style flags we had temporarily:**

```bash
-fcheck=all -finit-character=0 -finit-local-zero -ffpe-trap=invalid,zero,overflow -g -fbacktrace
```

**Resolution:**

Switch to production flags (remove runtime checking and aggressive traps). Then multi-thread runs pass.

**Why this matters:**

These debug flags are fine for single-thread diagnosis but can destabilize OpenMP execution. Always validate with `OMP_NUM_THREADS>1` before declaring success.

---

## Appendix B: Stack Overflow Investigation

### Symptom

Both QE `pw.exe` installed via MSYS2 pacman (package: `mingw-w64-ucrt-x86_64-quantum-espresso`) and QE `pw.exe` built locally following the MSYS2 recipe/PKGBUILD can appear to "hang" right after startup banner / printing program dimensions, even for a tiny smoke test (very small system, ~0.27s on other builds).

Example output stops after:
- "Serial multi-threaded version, running on 16 processor cores" (sometimes)
- "Reading input ..."
- "Current dimensions of program PWSCF are: ... lmaxx ..."

Then no progress.

### Investigation / Attempted Mitigations

First suspicion was OpenMP / library threads (FFTW/OpenBLAS/libgomp). Tried forcing single-thread and large OpenMP thread stacks:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_STACKSIZE=512M
export GOMP_STACKSIZE=512M
```

This did **NOT** fix the issue: still dies.

### Root Cause Proof

Running under strace shows the process terminates with Windows exception:

```
exception c00000fd
exited with status 0xc00000fd
```

Which corresponds to `STATUS_STACK_OVERFLOW` (stack overflow), i.e. it's not a hang.

Verified the PE header stack reserve is only 2 MB:

```bash
objdump -x "$(which pw.exe)" | grep -i -E "SizeOfStackReserve|SizeOfStackCommit"
```

Output:
```
SizeOfStackReserve 0000000000200000   (2MB)
SizeOfStackCommit  0000000000001000
```

This confirms the main executable stack is too small on Windows/MinGW for QE (Fortran code uses large stack frames / automatic arrays), so it overflows quickly, even on small cases.

### Final Conclusion

The MSYS2 pacman package as-is is **not reliable for QE on Windows for real workloads** because it links executables with an insufficient stack reserve (2MB).

The correct fix is to rebuild QE with a larger stack reserve at link time. Environment variables for OMP thread stack do **not** solve the main stack reserve problem.

### How to Fix

When building via PKGBUILD / recipe, add linker flag to increase stack reserve, e.g. 256MB:

```bash
export LDFLAGS+=" -Wl,--stack,268435456"
```

(128MB also may work: `134217728`, but recommend 256MB to be safe.)

After rebuild, verify:

```bash
objdump -x path/to/pw.exe | grep -i SizeOfStackReserve
```

Should show a much larger value (e.g., `0000000010000000` for 256MB).

### Quick Triage Checklist

If QE appears to "hang" early:

1. **Confirm you're using the intended pw.exe:**
   ```bash
   which pw.exe
   ```

2. **Capture strace to see the actual exception:**
   ```bash
   strace -f -o pw_strace.txt pw.exe -in input.in
   ```
   Look for `c00000fd` in the output.

3. **Check stack reserve with objdump** (requires binutils):
   ```bash
   pacman -S --needed mingw-w64-ucrt-x86_64-binutils
   objdump -x "$(which pw.exe)" | grep -i SizeOfStackReserve
   ```
   If it shows `0000000000200000` (2MB), that's the problem.

**Recommendation:** Do not rely on the stock MSYS2 QE package for production; rebuild with `--stack` flag.
