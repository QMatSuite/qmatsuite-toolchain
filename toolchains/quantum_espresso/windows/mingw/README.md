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
# Update package database
pacman -Syu

# Install required packages
pacman -S base-devel git wget patch \
  mingw-w64-x86_64-gcc \
  mingw-w64-x86_64-gcc-fortran \
  mingw-w64-x86_64-openblas \
  mingw-w64-x86_64-fftw \
  mingw-w64-x86_64-hdf5 \
  mingw-w64-x86_64-omp \
  mingw-w64-x86_64-pkg-config
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
export FFLAGS="-Wno-missing-include-dirs -w -std=legacy -fallow-argument-mismatch -fallow-invalid-boz"
export CFLAGS="-Wno-incompatible-pointer-types"
export QE_PREFIX="/c/Users/$(whoami)/qmatsuite-toolchain/toolchains/quantum_espresso/windows/mingw/upstream/qe/artifacts"
export FFT_CFLAGS="$(pkg-config --cflags fftw3)"
export FFT_LIBS="$(pkg-config --libs fftw3)"
# 可选：如果你想用 fftw threads
# export FFT_LIBS="$(pkg-config --libs fftw3 fftw3_threads)"

./configure \
  MPIF90=gfortran \
  BLAS_LIBS="-lopenblas" \
  LAPACK_LIBS="-lopenblas" \
  FFT_LIBS="$FFT_LIBS" \
  CPPFLAGS="$FFT_CFLAGS" \
  --with-scalapack=no \
  --enable-openmp \
  --prefix="$QE_PREFIX"
```

**Configuration options:**
- `F90=gfortran`: Use gfortran as Fortran compiler
- `MPIF90=gfortran`: Use gfortran for MPI (serial build)
- `BLAS_LIBS="-lopenblas"`: Link against OpenBLAS
- `LAPACK_LIBS="-lopenblas"`: Use OpenBLAS for LAPACK
- `FFT_LIBS="-lfftw3"`: Use FFTW3 for FFT
- `--with-scalapack=no`: Disable SCALAPACK (not available in MinGW)
- `--enable-openmp`: Enable OpenMP support

### Step 4: Build QE

Build the desired targets:

```bash
# Build all executables
# make -j$(nproc)
make -j1
# Or build specific targets (e.g., pw.x)
# should only use make pw -j1 due to competition problem
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
```

### Step 6: Stage Required DLLs

After building and installing, you need to copy the required DLLs from MinGW to your installation directory so that `pw.exe` can run. Use the staging script:

**From PowerShell:**
```powershell
cd toolchains\quantum_espresso\windows\mingw
.\scripts\stage_qe_dlls.ps1
```

**Or with custom paths:**
```powershell
.\scripts\stage_qe_dlls.ps1 -QeBinDir "C:\path\to\qe\bin" -MingwBinDir "C:\tools\msys64\mingw64\bin"
```

**Or manually copy DLLs in MSYS2 MinGW terminal:**
```bash
# Set paths
QE_BIN="/c/Users/$(whoami)/qmatsuite-toolchain/toolchains/quantum_espresso/windows/mingw/upstream/qe/artifacts/bin"
MINGW_BIN="/mingw64/bin"

# Copy required DLLs
cp "$MINGW_BIN/libgomp-1.dll" "$QE_BIN/"
cp "$MINGW_BIN/libgfortran"*.dll "$QE_BIN/"
cp "$MINGW_BIN/libopenblas"*.dll "$QE_BIN/"
cp "$MINGW_BIN/libfftw3"*.dll "$QE_BIN/"
cp "$MINGW_BIN/libgcc_s_seh-1.dll" "$QE_BIN/"
cp "$MINGW_BIN/libquadmath"*.dll "$QE_BIN/"
cp "$MINGW_BIN/libwinpthread-1.dll" "$QE_BIN/"
```

The staging script automatically:
- Detects your QE installation directory
- Finds your MSYS2/MinGW installation
- Copies all required DLLs
- Verifies that critical DLLs are present

The executables will be installed in `$PREFIX/bin/` with `.x` suffix (e.g., `pw.x`).

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

The following DLLs need to be available when running QE executables:

**Direct dependencies of pw.exe:**
- `libopenblas*.dll` - OpenBLAS runtime
- `libgfortran-*.dll` - GNU Fortran runtime (typically `libgfortran-5.dll`)
- `libgomp-1.dll` - OpenMP runtime (required if OpenMP enabled)
- `libgcc_s_seh-1.dll` - GCC runtime
- `libfftw3-3.dll` - FFTW3 runtime
- `libfftw3_omp-3.dll` - FFTW3 OpenMP support

**Indirect dependencies:**
- `libquadmath-*.dll` - Quadmath library (required by libgfortran)
- `libwinpthread-1.dll` - POSIX threads (required by libquadmath)

These are typically located in `/mingw64/bin/` (or `C:\tools\msys64\mingw64\bin` on Windows).

**Note:** Use the `stage_qe_dlls.ps1` script (see Step 6 above) to automatically copy all required DLLs to your installation directory.

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
