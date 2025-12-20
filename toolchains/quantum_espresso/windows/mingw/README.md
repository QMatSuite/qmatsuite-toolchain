# QE on Windows (MinGW)

Windows-specific build notes and scripts for Quantum ESPRESSO using MSYS2/MinGW-w64.

## Prerequisites

### 1. Install MSYS2

Install MSYS2 using Chocolatey:
```powershell
choco install msys2
```

Or download and install manually from: https://www.msys2.org/

### 2. Install Required Packages

Open **MSYS2 MinGW 64-bit terminal** (not the regular MSYS2 terminal) and run:

```bash


# Install required packages
pacman -Syu

# 基础工具（MSYS 里的，不分前缀）
pacman -S --needed base-devel git wget patch make pkgconf

# UCRT64 toolchain + libs（全部 ucrt 前缀）
pacman -S --needed \
  mingw-w64-ucrt-x86_64-toolchain \
  mingw-w64-ucrt-x86_64-gcc-fortran \
  mingw-w64-ucrt-x86_64-openblas \
  mingw-w64-ucrt-x86_64-fftw \
  mingw-w64-ucrt-x86_64-hdf5 \
  mingw-w64-ucrt-x86_64-pkgconf \
  mingw-w64-ucrt-x86_64-ntldd
```

## Building QE Locally

### Step 1: Download and Extract QE Source

In the MSYS2 MinGW 64-bit terminal:

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
if [ ! -d "qe" ]; then
    tar -xzf qe.tar.gz
    EXTRACTED=$(find . -maxdepth 1 -type d -name "q-e-${QE_VERSION}*" | head -1)
    if [ -n "$EXTRACTED" ]; then
        mv "$EXTRACTED" qe
    fi
fi

cd qe
```

### Step 2: Apply Patches (optional)

Apply the MinGW-specific patches from the repository:

```bash
# Navigate to QE source directory
cd /c/Users/$(whoami)/qmatsuite-toolchain/toolchains/quantum_espresso/windows/mingw/upstream/qe

# Set patch directory
PATCH_DIR="/c/Users/$(whoami)/qmatsuite-toolchain/toolchains/quantum_espresso/windows/mingw/patches"

# Apply patch 001: use-srand-instead-of-srandom
patch -Nbp1 -i "$PATCH_DIR/001-use-srand-instead-of-srandom.patch" || echo "Patch 001: may already be applied"

# Note: Patches 002-004 are for CMake builds and may not be needed for configure builds
# Apply them if you encounter issues:
# important!
patch -Nbp1 -i "$PATCH_DIR/002-fix-build-devxlib-on-mingw.patch"
# patch -Nbp1 -i "$PATCH_DIR/003-fix-install-c-libraries.patch"
patch -Nbp1 -i "$PATCH_DIR/004-fix-executable-suffix-on-mingw.patch"
```

### Step 3: Configure QE

Configure QE using the traditional `./configure` script:

```bash
# Set compiler flags to suppress warnings
export FFLAGS="-Wno-missing-include-dirs -w -std=legacy -fallow-argument-mismatch -fallow-invalid-boz -fopenmp -fmax-stack-var-size=0"
export CFLAGS="-O2 -Wno-incompatible-pointer-types"
export QE_PREFIX="/c/Users/$(whoami)/qmatsuite-toolchain/toolchains/quantum_espresso/windows/mingw/upstream/qe/artifacts"
#export FFT_CFLAGS="$(pkg-config --cflags fftw3)"
#export FFT_LIBS="$(pkg-config --libs fftw3)"
export LDFLAGS="-Wl,--stack,268435456 -fopenmp"
export FCFLAGS="-O2 -fopenmp -fmax-stack-var-size=0"

# 可选：如果你想用 fftw threads
# export FFT_LIBS="$(pkg-config --libs fftw3 fftw3_threads)"

./configure \
  MPIF90=gfortran \
  BLAS_LIBS="-lopenblas" \
  LAPACK_LIBS="-lopenblas" \
  FFT_LIBS="-lfftw3 -lfftw3_omp" \
  DFLAGS="-D_WIN32 -D__FFTW3" \
  --enable-openmp \
  --prefix="$QE_PREFIX"
```

**Configuration options:**
- The configuration should either detect external fftw3 automatically, or fall back to internal FFTW. Either way ok. Problem arises if user passes FFT_LIBS to configure (for example if one wants to pass -lfftw3_omp), then it will fail to decide, and neither -D__FFTW nor -D__FFTW3 flag will be set in DFLAGS, resulting in compilation error (cannot find some fftw include files). Therefore, the user should set -D__FFTW3 flag manually, signaling external FFTW3 usage. If the user mistakenly use -D__FFTW flag signaling internal FFTW usage instead, compilation will work, but there will be runtime problem of " Reading input from si100.in

     Current dimensions of program PWSCF are:
     Max number of different atomic species (ntypx) = 10
     Max number of k-points (npk) =  40000
     Max angular momentum in pseudopotentials (lmaxx) =  4
 end of file reached, tag  not found
 end of file reached, tag not found
At line 1081 of file xmltools.f90
Fortran runtime error: Read past ENDFILE record

Error termination. Backtrace:
#0  0xbd4e5c23 in ???
#1  0xbd461011 in ???
#2  0xbd31cafb in ???
#3  0xbd4626fd in ???
#4  0xbd4e33f1 in ???
#5  0xbd4e2dc9 in ???
#6  0xbd5005ac in ???
#7  0xbd4ffda3 in ???
#8  0xbd44d6e2 in ???
#9  0x6aa00c67 in ???
#10  0x6a9d2b1a in ???
#11  0x6a9d3ff6 in ???
#12  0x6a9beea7 in ???
#13  0x6a5dd889 in ???
#14  0x6a0d1784 in ???
#15  0x6a0d9321 in ???
#16  0x6a1abb1f in ???
#17  0x6a041614 in ???
#18  0x6a0416ba in ???
#19  0x6a0410c8 in __tmainCRTStartup
        at D:/W/B/src/mingw-w64/mingw-w64-crt/crt/crtexe.c:246
#20  0x6a041435 in mainCRTStartup
        at D:/W/B/src/mingw-w64/mingw-w64-crt/crt/crtexe.c:123
#21  0x4a77e8d6 in ???
#22  0x4aa6c53b in ???
#23  0xffffffff in ???
". This seemingly unrelated problem (reading xml problem) is actually related to missing/unmatched FFT lib.
- `F90=gfortran`: Use gfortran as Fortran compiler
- `MPIF90=gfortran`: Use gfortran for MPI (serial build)
- `BLAS_LIBS="-lopenblas"`: Link against OpenBLAS
- `LAPACK_LIBS="-lopenblas"`: Use OpenBLAS for LAPACK
- `--enable-openmp`: Enable OpenMP support

**Modify these lines in make.inc (might not be necessary)**

### Step 4: Build QE

Build the desired targets:

```bash
# Build all executables
# !!! should only use make pw -j1 due to competition problem
# make -j$(nproc)
make -j1
# Or build specific targets (e.g., pw.x)
# !!! should only use make pw -j1 due to competition problem
# make pw -j$(nproc)
make pw -j1
```

### Cleanly build QE (if things go wrong):
```bash
make veryclean 2>/dev/null || make clean
./configure ...（你原来的参数）...
make -j
```
### Step 5: Install (Optional)

Install to a prefix directory:

```bash
# Set installation prefix
# make install PREFIX="$QE_PREFIX"
make install 

export QE_BIN="$USERPROFILE/qmatsuite-toolchain/toolchains/quantum_espresso/windows/mingw/upstream/qe/artifacts/bin"
export PATH="$QE_BIN:$PATH"
```

### Step 6: Stage Executables and Required DLLs

After building and installing, you need to stage the executables and required DLLs to a clean distribution directory. Use the staging script:

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

1. **Finds executables**: Recursively searches for all `*.x` files in the artifacts directory (handles both `artifacts\bin` and `artifacts` root).

2. **Copies to distribution directory**: Copies all `*.x` files to `dist\win-mingw` in a flat layout (no subdirectories).

3. **Renames to .exe**: Renames all copied `*.x` files to `*.exe` in the distribution directory (e.g., `pw.x` → `pw.exe`).

4. **Copies required DLLs**: Copies required runtime DLLs from MSYS2 UCRT64 bin directory to the distribution directory.

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

## Running QE

### From MSYS2 MinGW Terminal

```bash
# Navigate to bin directory
cd /c/qe-install/bin

# Run QE
./pw.x < input.in > output.out
```

### From Windows Command Prompt

1. Navigate to the installation directory:
   ```cmd
   cd C:\qe-install\bin
   ```

2. Run QE:
   ```cmd
   pw.x < input.in > output.out
   ```

**Note:** Make sure all required DLLs are in the same directory as the executables, or add the MinGW bin directory to your PATH:
```cmd
set PATH=%PATH%;C:\tools\msys64\mingw64\bin
```

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

## Troubleshooting

### Build Errors

1. **Compiler not found**: Ensure you're using MSYS2 MinGW 64-bit terminal, not the regular MSYS2 terminal.

2. **Library not found**: Check that all required packages are installed:
   ```bash
   pacman -Qs mingw-w64-x86_64-openblas
   pacman -Qs mingw-w64-x86_64-fftw
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
# Ensure FFTW3 is installed in the same MSYS2 prefix
pacman -S mingw-w64-x86_64-fftw

# Use pkg-config to get correct flags
export FFT_CFLAGS="$(pkg-config --cflags fftw3)"
export FFT_LIBS="$(pkg-config --libs fftw3)"

# Reconfigure with explicit FFT flags
./configure \
  FFT_LIBS="$FFT_LIBS" \
  CPPFLAGS="$FFT_CFLAGS" \
  ... (other options)
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
export FFLAGS="-Wno-missing-include-dirs -w -std=legacy -fallow-argument-mismatch -fallow-invalid-boz"

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

1. **DLL not found**: Copy required DLLs to the executable directory or add MinGW bin to PATH.

2. **Path issues**: Use forward slashes in MSYS2 terminal, or use Windows paths with proper escaping.

## MSYS2/MinGW QE can crash with STACK_OVERFLOW (0xC00000FD) — fix is linker stack reserve

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

## Patches

The build process uses patches from the [MSYS2 MINGW-packages repository](https://github.com/msys2/MINGW-packages/tree/master/mingw-w64-quantum-espresso) to fix MinGW-specific build issues. These patches are stored locally in `toolchains/quantum_espresso/windows/mingw/patches/`.

See `toolchains/quantum_espresso/windows/mingw/patches/README.md` for details about each patch.

## CI Build (CMake)

The CI workflow (`.github/workflows/qe-windows-mingw.yml`) uses CMake instead of the configure script. For CMake-based builds, see the workflow file for the exact configuration options used.

## Migration from Intel oneAPI

This toolchain replaces the previous Intel oneAPI/ifort/MS-MPI based Windows builds. The MinGW-w64 toolchain provides:

- Open-source compiler (gfortran)
- No proprietary dependencies
- Better compatibility with open-source scientific software
- Easier redistribution (no Intel license restrictions)
