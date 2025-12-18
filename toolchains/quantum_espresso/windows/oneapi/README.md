# Quantum ESPRESSO Windows oneAPI Toolchain

This directory contains build scripts and documentation for building Quantum ESPRESSO on Windows using Intel oneAPI compilers (ifx/icx).

**Note on toolchain scope:** This document specifically covers the Intel oneAPI toolchain for Windows. There is a separate MinGW/MSYS2-based Windows toolchain for Quantum ESPRESSO in this repository under `toolchains/quantum_espresso/windows/mingw`. That MinGW path uses gfortran and MS-MPI, is already working, and we plan to compare its performance with the Intel oneAPI + MKL path. MinGW-specific details are documented separately in that toolchain's documentation.

## Quick Links

- **Build script**: `scripts/build_qe_win_oneapi.ps1`
- **Staging script**: `scripts/stage_qe_windows.ps1`
- **Refresh script**: `scripts/refresh_qe_source.ps1`
- **Patch script**: `scripts/apply_qe_patches.ps1`

## Why Git Clone is Required

**Important:** QE must be cloned from git, not downloaded as a zip archive. Here's why:

1. **Git Submodules**: QE uses git submodules (Wannier90, MBD, DeviceXlib) that are not included in zip archives
2. **CMake Requirements**: QE's CMake build system requires proper git metadata to initialize submodules correctly
3. **Version Tracking**: Git clones provide proper version information and commit history

GitHub/GitLab "Download ZIP" buttons do **not** include submodules, which will cause the build to fail with errors like:
- `CMake Error: C preprocessor cpp not found`
- `git fetch failed!`
- `Cannot find source file: wannier90/src/parameters.F90`

## Repository Structure

This repository uses a vendor/submodule approach:

```
.
├── scripts/                    # Our build scripts
│   ├── build_qe_win_oneapi.ps1 # Main build script
│   ├── refresh_qe_source.ps1   # QE source refresh script
│   ├── apply_qe_patches.ps1    # Windows patch application
│   └── stage_qe_windows.ps1     # Runtime DLL staging
├── upstream/                   # Vendor directory
│   └── qe/                    # QE source (cloned from git)
│       ├── PW/
│       ├── PHonon/
│       ├── external/          # Git submodules
│       └── ...
├── patches/                    # Windows-specific patches
└── dist/                       # Staged distributable binaries
```

## Prerequisites

1. **Intel oneAPI HPC Toolkit** (2024 or later)
   - Includes `ifx` (Fortran) and `icx` (C) compilers
   - Includes Intel MKL (BLAS/LAPACK)
   - Download from: https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html

2. **Visual Studio 2022 Build Tools** (or full Visual Studio)
   - Required for the MSVC linker (Intel compilers need it on Windows)
   - Install "Desktop development with C++" workload
   - Download from: https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2022

3. **CMake** (3.20 or later)
   - Download from: https://cmake.org/download/
   - Or install via Chocolatey: `choco install cmake`

4. **Ninja** (build system)
   - Install via Chocolatey: `choco install ninja`
   - Or download from: https://github.com/ninja-build/ninja/releases

5. **Git** (for cloning QE)
   - Download from: https://git-scm.com/download/win
   - Or install via Chocolatey: `choco install git`

## CI/CD: Automated Build and Staging

The GitHub Actions workflow (`.github/workflows/qe-windows-oneapi-mkl.yml`) automates the complete build and staging pipeline:

1. **Installs oneAPI** (with caching for performance)
2. **Builds QE** using the oneAPI toolchain
3. **Stages the distribution** using `stage_qe_windows.ps1` with `-Strict` mode
4. **Uploads the staged distribution** as a GitHub Actions artifact (`qe-windows-oneapi-dist`)

The staged artifact contains:
- All QE executables (`.exe` files)
- Required runtime DLLs (Intel Fortran, OpenMP, MKL, MSVC runtime)
- Test resources (`scf-cg.in`, `Si.pz-vbc.UPF`)
- Dependency closure log (`deps-closure.txt`)

The artifact is retained for 7 days and can be downloaded from the workflow run page.

## CI/CD Installation: oneAPI with MKL on Windows

This section documents our experience installing Intel oneAPI HPC Toolkit (with MKL) in GitHub Actions Windows CI environments. The official oneAPI CI samples repository ([oneapi-src/oneapi-ci](https://github.com/oneapi-src/oneapi-ci)) and [documentation](https://oneapi-src.github.io/oneapi-ci/#intelr-oneapi-hpc-toolkit) provide component listings and installation methods, but we encountered several issues that required a specific approach.

### What We Needed

For Quantum ESPRESSO builds, we require:
- **Intel Fortran compiler** (`ifx`/`ifort`) - Component ID: `intel.oneapi.win.ifort-compiler`
- **Intel C++ compiler** (`icx`) - Component ID: `intel.oneapi.win.cpp-dpcpp-common`
- **Intel MKL** (BLAS/LAPACK) - Component ID: `intel.oneapi.win.mkl.devel`

The MKL component is critical; without it, CMake's `FindLAPACK` cannot locate BLAS/LAPACK libraries, causing build failures.

### Installation Methods We Tried (and Why They Failed)

#### 1. **winget (Windows Package Manager)** ❌

**Attempt:**
```powershell
winget install Intel.oneAPI.HPCKit
```

**Why it failed:**
- winget installs a minimal oneAPI installation that **does not include MKL**
- The installed components are stripped down and missing critical libraries
- CMake configuration fails with "LAPACK not found" errors
- No way to specify which components to install via winget

#### 2. **Chocolatey** ❌

**Attempt:**
```powershell
choco install intel-oneapi-hpckit
```

**Why it failed:**
- Similar to winget, Chocolatey packages often install minimal oneAPI distributions
- MKL components are frequently missing or incomplete
- Component selection is limited or non-existent

#### 3. **Direct Installer with `--components` Flag** ❌

**Attempt:**
```powershell
.\intel-oneapi-hpc-toolkit-2025.3.0.378_offline.exe --components intel.oneapi.win.cpp-dpcpp-common:intel.oneapi.win.ifort-compiler:intel.oneapi.win.mkl.devel
```

**Why it failed:**
- On Windows CI runners, the installer **hangs indefinitely** when using `--components` flag directly
- The installer appears to start but never completes (waits forever)
- This is a known issue with the oneAPI installer on Windows in automated environments
- No timeout mechanism; the CI job hangs until manually cancelled

#### 4. **Web Installer** ❌

**Why it failed:**
- Requires interactive GUI, not suitable for CI
- Cannot be automated in headless environments

### The Working Solution: Offline Installer with Config File ✅

Based on the [oneAPI CI samples repository](https://github.com/oneapi-src/oneapi-ci) and community discussions, the solution is to:

1. **Extract the offline installer** (silent extraction)
2. **Create an installation config file** (`.ini` format) specifying components
3. **Run the bootstrapper** with the config file (not the main installer)

**Why this works:**
- The bootstrapper (`bootstrapper.exe`) handles component selection reliably
- The config file approach avoids the hanging issue with `--components` flag
- Silent installation works in CI environments
- All specified components (including MKL) are installed correctly

**Implementation:**

```powershell
# 1. Download offline installer
$url = "https://registrationcenter-download.intel.com/akdlm/IRC_NAS/3a871580-f839-46ed-aeae-685084127279/intel-oneapi-hpc-toolkit-2025.3.0.378_offline.exe"
$installerExe = "oneapi_hpc_installer.exe"
Invoke-WebRequest -Uri $url -OutFile $installerExe

# 2. Extract installer silently
$extractDir = "oneapi_extracted"
Start-Process -FilePath $installerExe -ArgumentList "-s", "-x", "-f", $extractDir -Wait

# 3. Create config file
$configFile = "oneapi_install_config.ini"
@"
eula=accept
components=intel.oneapi.win.cpp-dpcpp-common:intel.oneapi.win.ifort-compiler:intel.oneapi.win.mkl.devel
p="NEED_VS2017_INTEGRATION=0 NEED_VS2019_INTEGRATION=0 NEED_VS2022_INTEGRATION=0"
"@ | Out-File -FilePath $configFile

# 4. Run bootstrapper with config file
$bootstrapper = Join-Path $extractDir "bootstrapper.exe"
$installLogDir = "install_logs"
Start-Process -FilePath $bootstrapper -ArgumentList "--silent", "--eula", "accept", "--components", "intel.oneapi.win.cpp-dpcpp-common:intel.oneapi.win.ifort-compiler:intel.oneapi.win.mkl.devel", "--log-dir", $installLogDir -Wait
```

**Key points:**
- Component IDs must match exactly those listed in the [oneAPI CI documentation](https://oneapi-src.github.io/oneapi-ci/#intelr-oneapi-hpc-toolkit)
- The `p="..."` parameter disables Visual Studio integration prompts (not needed in CI)
- The bootstrapper approach is the only reliable method for automated Windows CI installation

**Reference implementation:** See `.github/workflows/qe-windows-oneapi-mkl.yml` in this repository for a complete GitHub Actions workflow that implements this installation method with caching, builds QE, stages the distribution, and uploads it as an artifact.

### Caching for CI Performance

Since oneAPI installation takes 20+ minutes, we cache both:
1. **The installer download** (several GB) - cached in `.github/cache/oneapi-installer/`
2. **The installed components** - cached in `C:\Program Files (x86)\Intel\oneAPI\`

Cache keys include version numbers so updates trigger fresh installations:
```yaml
- uses: actions/cache@v4
  with:
    path: |
      C:\Program Files (x86)\Intel\oneAPI\compiler
      C:\Program Files (x86)\Intel\oneAPI\mkl
      C:\Program Files (x86)\Intel\oneAPI\setvars.bat
    key: oneapi-install-2025.3.0-${{ runner.os }}-cpp-dpcpp-common-ifort-compiler-mkl-devel
```

### Component ID Reference

From the [official oneAPI CI documentation](https://oneapi-src.github.io/oneapi-ci/#intelr-oneapi-hpc-toolkit), the correct component IDs for Windows HPC Toolkit 2025.3.0 are:

| Component | ID | Provides |
|-----------|-----|----------|
| C++ Compiler | `intel.oneapi.win.cpp-dpcpp-common` | `icx`, `icpx`, `dpcpp` |
| Fortran Compiler | `intel.oneapi.win.ifort-compiler` | `ifx`, `ifort` |
| MKL Development | `intel.oneapi.win.mkl.devel` | MKL headers, libraries, CMake configs |

**Important:** Use `intel.oneapi.win.mkl.devel` (not `intel.oneapi.win.mkl`) to get CMake package configs that `FindLAPACK` can locate.

### Verifying Installation

After installation, verify MKL is accessible:
```powershell
# Initialize oneAPI environment
& "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64

# Check MKLROOT is set
echo $env:MKLROOT

# Verify CMake can find MKL
cmake -DCMAKE_PREFIX_PATH="$env:MKLROOT\lib\cmake\mkl" -P FindMKL.cmake
```

The installation should be at: `C:\Program Files (x86)\Intel\oneAPI\`

## Quick Start (end-to-end workflow)

Order matters: refresh (clone) → apply patches → build → stage.

### Step 0: Install prerequisites
- Intel oneAPI HPC Toolkit (2024+). Includes ifx/icx, MKL, setvars.bat.
- Visual Studio 2022 (or 2019) Build Tools with "Desktop development with C++". oneAPI supports VS2019/VS2022; older versions do not work.
- Chocolatey (recommended) plus packages:
  ```powershell
  choco install git
  choco install cmake
  choco install ninja
  ```
  (If you install CMake/Ninja manually, ensure they're on PATH.)

### Step 1: Refresh (clone QE from git)
```powershell
.\scripts\refresh_qe_source.ps1 -NoPatch
```
Notes:
- Cleans prior QE sources/build artifacts, preserves repo files.
- Clones QE from GitLab with submodules into `upstream/qe`.
- To target a tag: `.\scripts\refresh_qe_source.ps1 -QeVersion qe-7.5`
- Dry run: `.\scripts\refresh_qe_source.ps1 -DryRun`

### Step 2: Apply Windows patches
```powershell
.\scripts\apply_qe_patches.ps1
```
Notes:
- Applies `qe-win-cmake-generation` and `qe-win-c-portability` patches to `upstream/qe`.
- Replaces `external/devxlib/src/timer.c` with Windows-compatible version from `patches/devxlib-timer.c`.
- Idempotent: skips patches if already applied; file replacement always overwrites.

### Step 3: Build QE (oneAPI + Ninja)
```powershell
.\scripts\build_qe_win_oneapi.ps1
```
Options:
- Serial build: `-NoMpi`
- Custom source dir: `-QeSourceDir "path/to/qe"`

### Step 4: Stage runtime-ready binaries
```powershell
.\scripts\stage_qe_windows.ps1
```
Notes:
- Outputs to `dist/win-oneapi/`
- Implements three-layer staging pipeline:
  1. **Recursive dependency closure** - Uses `dumpbin /DEPENDENTS` to discover all transitive DLL dependencies (auto-discovers `dumpbin.exe` from VS2022/2019)
  2. **Must-have DLL patterns** - Ensures critical runtime DLLs (Fortran, OpenMP, MKL, MSVC) are included
  3. **Clean-room smoke test** - Runs `pw.exe -i scf-cg.in` and verifies "JOB DONE" in output
- Filter executables: `-Only pw.x,ph.x`
- Clean staged dir: `-Clean`
- Strict mode (fail if smoke test doesn't pass): `-Strict`
- Generates `deps-closure.txt` with all resolved DLLs
- Copies test resources from `toolchains/quantum_espresso/tests/resources/` for smoke test

## Build Process

The build script performs these steps:

1. **Verifies QE source and submodules** - Ensures QE is properly cloned with submodules initialized
2. **Discovers Intel oneAPI** - Finds `setvars.bat` to initialize compiler environment
3. **Checks Visual Studio** - Verifies MSVC linker is available
4. **Locates CMake and Ninja** - Finds build tools
5. **Configures CMake** - Sets up build with Intel compilers and MKL
6. **Builds QE** - Compiles `pw.x`, `ph.x`, `pp.x` executables

## Build Output

After a successful build, executables are located in:
```
upstream/qe/build-win-oneapi/bin/
```

You should find:
- `pw.x` - Plane-wave self-consistent field calculations
- `ph.x` - Phonon calculations
- `pp.x` - Post-processing utilities
- And other QE executables

## Windows patch workflow

QE upstream source is not committed; Windows fixes live as patches under `patches/`. This repository carries small Windows-focused patches for the upstream QE sources.

### Patch List

- **`qe-win-cmake-generation.patch`**  
  Makes QE CMake generation Windows-friendly by avoiding shell pipelines and using CMake scripting for git info and Fortran includes. Includes:
  - CMake/build helpers
  - Fortran include preprocessing without shell redirection
  - Git revision header generation without sed/pipes

- **`qe-win-c-portability.patch`**  
  Adds Windows stubs and compatibility fixes for socket and math code. Includes:
  - Windows socket implementation (`Modules/sockets.c`)
  - Math constant definitions (`XClib/beefun.c`)

- **`devxlib-timer.c`** (file replacement, not a patch)  
  Complete file replacement for `external/devxlib/src/timer.c` with Windows-specific implementation using `QueryPerformanceCounter` instead of POSIX `gettimeofday()`. This replaces the original `qe-win-devxlib-timer.patch` approach (see "The `timer.c` Patching Struggle" below for details).

### Timer.c File Reference

The `patches/` directory contains reference files for the `timer.c` modification:

- **`patches/reference/timer.c.original`** - The original Unix version from QE upstream that uses `gettimeofday()` (POSIX-only, doesn't work on Windows)
- **`patches/devxlib-timer.c`** - The Windows-compatible replacement that uses `QueryPerformanceCounter` on Windows and falls back to `gettimeofday()` on Unix/Linux

The patching script (`apply_qe_patches.ps1`) copies `patches/devxlib-timer.c` directly to `upstream/qe/external/devxlib/src/timer.c`, replacing the original file. The original is kept in `patches/reference/` for reference and comparison.

### How to Apply Patches

**Automatic (recommended):**
```powershell
.\scripts\refresh_qe_source.ps1
```
The refresh script applies patches automatically (opt-out with `-NoPatch`).

**Manual apply (if you clone QE yourself):**
```powershell
.\scripts\apply_qe_patches.ps1 -QeDir upstream/qe -PatchDir patches
```

The patch script is **idempotent**: it skips patches already applied and stops if a patch cannot be applied cleanly.

**Note on CI caching:** In CI environments (like GitHub Actions), the QE source directory is often cached. If the cached QE source already has patches applied, the script will detect this and skip re-applying them. The script checks for already-applied state first (before attempting to apply) to avoid confusing error messages. This is expected behavior and not an error.

### Workflow Reminder

The typical workflow is:
1. **Clone QE source:** `scripts\refresh_qe_source.ps1` only cleans and clones QE (use `-NoPatch` to skip patching)
2. **Apply patches:** Run `scripts\apply_qe_patches.ps1` to apply the Windows patches above
3. **Build:** Use `scripts\build_qe_win_oneapi.ps1` (add `-NoMpi` for serial builds)

### Reviewing Applied Patches

To see what's been patched:
```powershell
cd upstream/qe
git status
git diff --stat
```

### The `timer.c` Patching Struggle

The `external/devxlib/src/timer.c` file required Windows-specific changes to use `QueryPerformanceCounter` instead of Unix `gettimeofday()`. Our initial approach was to create a patch file (`qe-win-devxlib-timer.patch`), but this led to significant issues:

#### Local Development Issues

1. **Patch path problems:** The patch file initially had incorrect paths relative to the QE root directory
2. **Line ending issues:** Windows vs Unix line endings (`CRLF` vs `LF`) caused patch application failures
3. **Multiple rebuild attempts:** We rebuilt the patch multiple times with:
   - Correct paths: `external/devxlib/src/timer.c` (relative to QE root)
   - Normalized Unix line endings
   - Verified with `git apply --check` from QE root

**Local solution:** After fixing paths and line endings, the patch passed `git apply --check` and appeared to work locally.

#### CI Failure

Despite working locally, the patch **failed in GitHub Actions CI** with:
```
error: patch failed: external/devxlib/src/timer.c
error: external/devxlib/src/timer.c: patch does not apply
error: corrupted patch
```

**Root cause:** Even with normalized line endings and correct paths, the patch format itself was fragile. The CI environment's git configuration, line ending handling, or patch application mechanism differed from local development, causing the patch to be rejected as "corrupted."

#### Final Solution: Complete File Replacement

Rather than continuing to fight with patch format issues, we switched to a **complete file replacement** approach:

1. **Removed the patch file:** Deleted `patches/qe-win-devxlib-timer.patch`
2. **Preserved original for reference:** Moved the original Unix version to `patches/reference/timer.c.original`
3. **Created replacement file:** Added `patches/devxlib-timer.c` containing the complete Windows-compatible implementation (uses `QueryPerformanceCounter` on Windows, `gettimeofday()` on Unix)
4. **Updated patch script:** Modified `apply_qe_patches.ps1` to copy the replacement file directly:
   ```powershell
   Copy-Item -Path "$PatchDir\devxlib-timer.c" -Destination "$QeDir\external\devxlib\src\timer.c" -Force
   ```

**Why this works:**
- No patch format to parse or apply
- No line ending issues (PowerShell handles file copying correctly)
- Works identically in local development and CI
- Simpler and more reliable than patching

**Lesson learned:** For complete file replacements (especially small files), direct file copying is more reliable than patch files in CI environments where patch application can be sensitive to git configuration, line endings, and environment differences.

## Windows: staging distributable binaries

After a build, stage the executables and required Intel runtimes into a distributable folder:
```powershell
.\scripts\stage_qe_windows.ps1
```
- Output goes to `dist/win-oneapi/`
- Copy only specific executables (names with or without `.exe`):
```powershell
.\scripts\stage_qe_windows.ps1 -Only pw.x,ph.x
```
- Clean the staged folder before re-staging:
```powershell
.\scripts\stage_qe_windows.ps1 -Clean
```
- Use strict mode (fail if smoke test doesn't pass):
```powershell
.\scripts\stage_qe_windows.ps1 -Strict
```

The script implements a three-layer staging pipeline:

1. **Recursive dependency closure** - Uses `dumpbin /DEPENDENTS` to discover all transitive DLL dependencies (if `dumpbin.exe` is available)
2. **Must-have DLL patterns** - Ensures critical runtime DLLs (Fortran, OpenMP, MKL, MSVC) are included even if not detected by dumpbin
3. **Clean-room smoke test** - Runs a minimal QE calculation (`pw.exe -i scf-cg.in`) to verify the staged distribution works correctly

The script copies QE `.exe` files plus Intel OpenMP, MKL, and Fortran runtime DLLs from your oneAPI install so the binaries run on machines without the full toolchain.

**Output files:**
- `deps-closure.txt` - Tab-separated list of all resolved DLLs (`<basename>\t<full_path>`)
- Test resources (`scf-cg.in`, `Si.pz-vbc.UPF`) copied from `toolchains/quantum_espresso/tests/resources/`

**Note:** In CI environments (GitHub Actions), the staging step runs automatically after the build and uploads the staged distribution as an artifact. See the "CI/CD: Automated Build and Staging" section above for details.

To force dumpbin discovery, ensure `VCToolsInstallDir` is set, e.g.:
```powershell
$env:VCToolsInstallDir="C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\"
.\scripts\stage_qe_windows.ps1
```

## CMake / Windows workarounds applied

These are the key fixes that made the Windows/oneAPI CMake flow reliable:

- **Fortran include preprocessing without shells**  
  Added `qe_prepare_fortran_include` in `cmake/qeHelpers.cmake` to generate `.fh/.inc` files without `>` redirection or Unix pipes. It:
  - Detects `#`-directives; copies verbatim when none.  
  - Otherwise preprocesses with the Fortran compiler (`ifx/ifort: -fpp -E -P`, GNU: `-cpp -E -P`, NVHPC: `-Mpreprocess -E`) via a CMake `-P` script using `OUTPUT_FILE`, so Ninja/MSBuild work on Windows.  
  - Reuses QE global `-D` definitions.  
  LAXlib now calls this helper for all its headers (no maintain/allow lists).

- **C preprocessor choice**  
  QE expects `QE_CPP`; we pass `icx` explicitly to avoid missing `cpp` on Windows.

- **Git revision header without sed/pipes**  
  `cmake/GitInfo.cmake` now uses pure `execute_process` + `file(WRITE)` to create `git-rev.h`, with fallbacks to `UNKNOWN` when git is unavailable. No `sed`, `cut`, or shell redirection.

- **Submodules and internal deps**  
  Always clone QE via `scripts/refresh_qe_source.ps1` (git + submodules). In CMake configure we force internal FFTW (`-DQE_FFTW_VENDOR=Internal`) and internal Wannier90/MBD/DeviceXlib to avoid CMake-driven submodule fetches on Windows.

- **Headers and preprocessing order**  
  LAXlib headers are generated before compilation in `build-win-oneapi/LAXlib/include`; the helper ensures the output directory exists and uses Fortran preprocessing where needed.

- **Windows source fixes (already in tree)**  
  - `Modules/sockets.c` uses Winsock headers on Windows.  
  - `external/devxlib/src/timer.c` uses `QueryPerformanceCounter` for high-resolution timing on Windows.  
  - `XClib/beefun.c` defines `M_PI` and includes `<math.h>`.

- **Clean configure/build re-entry**  
  Build script removes stale `CMakeCache.txt` before reconfigure to avoid generator mismatches.

If you re-run CMake manually, keep these options consistent (matching `build_qe_win_oneapi.ps1`):
```
-G Ninja
-DCMAKE_Fortran_COMPILER=ifx
-DCMAKE_C_COMPILER=icx
-DQE_CPP=".../icx.exe"
-DQE_ENABLE_MPI=OFF|ON
-DQE_LAPACK_INTERNAL=OFF
-DQE_WANNIER90_INTERNAL=ON -DQE_MBD_INTERNAL=ON -DQE_DEVICEXLIB_INTERNAL=ON
-DQE_FFTW_VENDOR=Internal
```

## Troubleshooting

### "QE source directory not found"

**Solution:** Run `.\scripts\refresh_qe_source.ps1` first to clone QE.

### "Required submodules are missing"

**Cause:** QE was downloaded as a zip archive instead of cloned from git.

**Solution:**
1. Delete the QE source directory
2. Run `.\scripts\refresh_qe_source.ps1` to clone from git
3. Or manually: `cd upstream/qe && git submodule update --init --recursive`

### "CMake Error: C preprocessor cpp not found"

**Cause:** Submodules are not initialized.

**Solution:** Ensure QE is cloned from git with submodules. Run the refresh script.

### "icx: error #10037: could not find 'link'"

**Cause:** Visual Studio Build Tools not found or not initialized.

**Solution:**
1. Install Visual Studio 2022 Build Tools with "Desktop development with C++" workload
2. Ensure the script can find VS at standard locations

### "Ninja not found"

**Solution:** Install Ninja:
```powershell
choco install ninja
```

Or download from https://github.com/ninja-build/ninja/releases and add to PATH.

### Build fails with linker errors

**Cause:** Visual Studio environment not properly initialized.

**Solution:** The script should automatically initialize VS. If it fails, manually run:
```cmd
"C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\Tools\VsDevCmd.bat" -arch=x64
```

Then run the build script again.

### "LAPACK not found" or "MKL not found" in CMake

**Cause:** Intel MKL was not installed or CMake cannot locate it.

**Solutions:**
1. **Verify MKL is installed:**
   ```powershell
   Test-Path "C:\Program Files (x86)\Intel\oneAPI\mkl\latest"
   ```

2. **Check component installation:** If using CI/automated installation, ensure you used the bootstrapper method with config file (see "CI/CD Installation" section above). winget/Chocolatey installations often exclude MKL.

3. **Set MKLROOT and CMAKE_PREFIX_PATH:**
   ```powershell
   $env:MKLROOT = "C:\Program Files (x86)\Intel\oneAPI\mkl\latest"
   $env:CMAKE_PREFIX_PATH = "$env:MKLROOT\lib\cmake\mkl;$env:CMAKE_PREFIX_PATH"
   ```

4. **Verify CMake can find MKL:**
   ```powershell
   cmake -DCMAKE_PREFIX_PATH="$env:MKLROOT\lib\cmake\mkl" -P FindMKL.cmake
   ```

**Note on `-DBLA_VENDOR` flag:** The build script does not use `-DBLA_VENDOR=Intel10_64lp_seq` to force FindLAPACK to use Intel MKL. This flag was originally introduced because we thought MKL was installed but FindLAPACK couldn't locate it. However, we later discovered that MKL wasn't actually installed at that time. Once MKL is properly installed (with the `intel.oneapi.win.mkl.devel` component) and `CMAKE_PREFIX_PATH` is set to include MKL's CMake package directory, FindLAPACK automatically detects and uses MKL without needing the vendor flag. We removed the flag to keep the build configuration more generalized and compatible with other BLAS/LAPACK implementations.

### Why We Use External MKL Instead of Internal LAPACK

This toolchain uses external Intel MKL for BLAS/LAPACK (`-DQE_LAPACK_INTERNAL=OFF`) rather than QE's internal LAPACK implementation. This was a deliberate design decision based on our experience:

In theory, `QE_LAPACK_INTERNAL=ON` should work on Windows with oneAPI compilers. However, local experiments revealed configure failures and compatibility issues between QE's internal LAPACK and the Intel oneAPI toolchain. Based on prior experience, mixing GNU-compiled LAPACK (which QE's internal LAPACK typically uses) with Intel compilers on Windows tends to be fragile and can lead to runtime issues.

Therefore, we decided to standardize on external Intel MKL for BLAS/LAPACK on this oneAPI path and not invest more time into the internal LAPACK route. This approach provides:
- Better performance (MKL is highly optimized for Intel architectures)
- Proven compatibility with Intel compilers
- Simplified build configuration (no need to compile LAPACK from source)
- Consistent behavior across local development and CI environments

### Runtime crash: forrtl severe (170) stack overflow (pw.exe)

#### Symptom

Build succeeds, but `pw.exe` crashes immediately on simple input with `forrtl: severe (170): stack overflow`. The crash happens right after PWSCF prints the "Current dimensions…" banner.

Example crash output:

```
Reading input from .\scf-cg.in

     Current dimensions of program PWSCF are:

     Max number of different atomic species (ntypx) = 10        

     Max number of k-points (npk) =  40000

     Max angular momentum in pseudopotentials (lmaxx) =  4      

forrtl: severe (170): Program Exception - stack overflow

Image              PC                Routine            Line        Source

pw.exe             00007FF7CEF20D07  Unknown               Unknown  Unknown

pw.exe             00007FF7CE60C006  Unknown               Unknown  Unknown

pw.exe             00007FF7CE526C0B  Unknown               Unknown  Unknown

pw.exe             00007FF7CE5250D2  Unknown               Unknown  Unknown

pw.exe             00007FF7CE521190  Unknown               Unknown  Unknown

pw.exe             00007FF7CEF17B2B  Unknown               Unknown  Unknown

pw.exe             00007FF7CEF20C30  Unknown               Unknown  Unknown

KERNEL32.DLL       00007FFA4A77E8D7  Unknown               Unknown  Unknown

ntdll.dll          00007FFA4AA6C53C  Unknown               Unknown  Unknown
```

#### Root Cause

This is a **stack overflow**, not floating-point overflow. Windows default thread stack is small (typically 1MB), and QE/Fortran can allocate large automatic arrays and temporary arrays on the stack. This triggers the crash early in runtime when these arrays are allocated.

#### Fix

The fix is implemented in `scripts/build_qe_win_oneapi.ps1` and includes:

1. **Use Intel Fortran `-heap-arrays` flag** - Forces large automatic arrays and temporaries to be allocated on the heap instead of the stack
2. **Increase Windows executable stack reserve** - Set linker flag `/STACK:134217728` (128MB) to provide more stack space
3. **Add `-traceback` in Debug builds** - Improves diagnostic backtraces for debugging

The relevant CMake configure flags are:

```cmake
-DCMAKE_Fortran_FLAGS_RELEASE="-heap-arrays"
-DCMAKE_Fortran_FLAGS_DEBUG="-heap-arrays -traceback"
-DCMAKE_EXE_LINKER_FLAGS="/STACK:134217728"
```

**Why 128MB:** We confirmed 256MB also works, and it generally has no practical downside because stack reserve is mostly virtual address space. We chose 128MB as a more conservative, still-safe default for distribution; it's large enough for QE while avoiding an unnecessarily large reserve. The `-heap-arrays` flag is the primary fix; the stack reserve is a defensive safeguard.

**Note:** This change is intentionally minimal and does not alter MKL/BLAS/LAPACK discovery behavior. The build script still uses the same MKL detection logic as before.

#### Verification

After applying these flags, `pw.exe` should pass the minimal `scf-cg.in` test without the immediate stack overflow crash.

### Runtime error: Missing MKL DLLs (mkl_avx2.2.dll, mkl_def.2.dll)

#### Symptom

After staging executables, running `pw.exe` fails with:

```
INTEL oneMKL ERROR: The specified module could not be found. mkl_avx2.2.dll.
Intel oneMKL FATAL ERROR: Cannot load mkl_avx2.2.dll or mkl_def.2.dll
```

#### Root Cause

The staging script (`scripts/stage_qe_windows.ps1`) was not copying all required MKL interface DLLs. While it copied core MKL DLLs (`mkl_core*.dll`, `mkl_intel_lp64*.dll`, etc.), it was missing the AVX2-optimized interface DLLs (`mkl_avx2*.dll`) and the default interface DLL (`mkl_def*.dll`) that MKL dynamically loads at runtime.

#### Fix

The fix is implemented in `scripts/stage_qe_windows.ps1` by adding the missing DLL patterns to the MKL pattern list:

**Location 1** (line 406 - main MKL pattern copy):
```powershell
Copy-DllPatterns -Patterns @("mkl_core*.dll","mkl_intel_lp64*.dll","mkl_intel_thread*.dll","mkl_rt*.dll","mkl_avx2*.dll","mkl_def*.dll") ...
```

**Location 2** (line 415 - supplemental pattern copy):
```powershell
Copy-DllPatterns -Patterns @("mkl_intel_lp64*.dll","mkl_avx2*.dll","mkl_def*.dll","libintlc*.dll") ...
```

These patterns will match:
- `mkl_avx2.2.dll` (the specific DLL from the error)
- `mkl_def.2.dll` (the specific DLL from the error)
- Any future versions (e.g., `mkl_avx2.3.dll`, `mkl_def.3.dll`)

The `Get-SearchRoots` function already searches in MKL directories (`mkl\latest\bin`, `mkl\latest\redist\intel64`, `mkl\latest\redist\intel64\mkl`), so these DLLs will be found and copied to the dist directory.

#### Verification

After re-running the staging script, the dist directory should contain `mkl_avx2.2.dll` and `mkl_def.2.dll`, and `pw.exe` should run without the MKL DLL loading error.

## Advanced Usage

### Building Specific Targets

After CMake configuration, you can build specific targets:
```powershell
cd upstream/qe/build-win-oneapi
cmake --build . --target pw --parallel
```

### Custom CMake Options

Edit `scripts/build_qe_win_oneapi.ps1` to add custom CMake options, or configure manually:
```powershell
cd upstream/qe
cmake -S . -B build-custom -G Ninja -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_C_COMPILER=icx ...
```

### Cleaning Build Artifacts

To clean build artifacts while keeping QE source:
```powershell
Remove-Item -Recurse -Force upstream/qe/build-win-oneapi
```

To completely refresh QE source (clean and re-clone):
```powershell
.\scripts\refresh_qe_source.ps1
```

## Script Reference

### `scripts/refresh_qe_source.ps1`

Refreshes QE source by cleaning old sources and cloning fresh from git.

**Parameters:**
- `-DryRun` - Show what would be deleted without actually doing it
- `-QeVersion <tag>` - Clone specific version (e.g., "qe-7.5"). Default: newest release
- `-NoPatch` - Skip automatic patch application

**What it preserves:**
- `scripts/` directory
- `.github/` and other CI configs
- `docs/` directory
- Project config files (`.gitignore`, `.editorconfig`, `LICENSE`, etc.)

**What it removes:**
- QE source directories (PW/, PHonon/, Modules/, etc.)
- `external/` submodules
- Build artifacts (build/, CMakeCache.txt, etc.)

### `scripts/build_qe_win_oneapi.ps1`

Builds QE using Intel oneAPI compilers.

**Parameters:**
- `-NoMpi` - Build serial QE without MPI support
- `-QeSourceDir <path>` - Path to QE source (default: "upstream/qe")

**Environment Variables:**
- `ONEAPI_ROOT` - Override oneAPI installation path (default: `C:\Program Files (x86)\Intel\oneAPI`)

### `scripts/stage_qe_windows.ps1`

Stages QE executables and required runtime DLLs into a distributable directory. Implements a three-layer staging pipeline:

1. **Recursive dependency closure** - Uses `dumpbin /DEPENDENTS` to discover all transitive DLL dependencies
2. **Must-have DLL patterns** - Ensures critical runtime DLLs (Fortran, OpenMP, MKL, MSVC) are included even if not detected by dumpbin
3. **Clean-room smoke test** - Runs a minimal QE calculation to verify the staged distribution works correctly

**Parameters:**
- `-BuildDir <path>` - Path to QE build directory (default: "upstream/qe/build-win-oneapi")
- `-DistDir <path>` - Target dist directory (default: "dist/win-oneapi")
- `-Only <list>` - List of executable basenames to stage (e.g., "pw.x,ph.x")
- `-Clean` - Delete dist directory before staging
- `-VerboseDeps` - Print parsed dependencies for each executable/DLL
- `-WhatIf` - Dry run: print actions without copying/deleting
- `-ProbeExe <name>` - Executable name to probe after staging (default: pw.x.exe if present)
- `-Strict` - Fail if smoke test exit code != 0 OR "JOB DONE" not detected in output

**Output:**
- Staged executables and DLLs in `$DistDir`
- `deps-closure.txt` - Tab-separated list of all resolved DLLs (`<basename>\t<full_path>`)
- Test resources copied from `toolchains/quantum_espresso/tests/resources/` (for smoke test)

**Smoke Test:**
The script runs `pw.exe -i scf-cg.in` in the staged directory and checks for "JOB DONE" in the output. With `-Strict` mode, the script fails if the test doesn't pass, ensuring the staged distribution is functional.

## Additional Resources

- **QE Official Documentation:** https://www.quantum-espresso.org/
- **Intel oneAPI Documentation:** https://www.intel.com/content/www/us/en/developer/tools/oneapi/documentation.html
- **QE GitHub Repository:** https://github.com/QEF/q-e
- **Intel oneAPI CI Samples:** https://github.com/oneapi-src/oneapi-ci
  - Sample CI configurations for installing oneAPI in various CI systems
  - Component listings and installation methods
- **oneAPI CI Component Documentation:** https://oneapi-src.github.io/oneapi-ci/#intelr-oneapi-hpc-toolkit
  - Complete list of component IDs for Windows/Linux/macOS
  - Version information and component dependencies
