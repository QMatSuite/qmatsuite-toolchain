# Windows Build Guide for Quantum ESPRESSO

This guide explains how to build Quantum ESPRESSO (QE) on Windows using Intel oneAPI compilers.

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
│   └── refresh_qe_source.ps1   # QE source refresh script
├── upstream/                   # Vendor directory
│   └── qe/                    # QE source (cloned from git)
│       ├── PW/
│       ├── PHonon/
│       ├── external/          # Git submodules
│       └── ...
└── docs/                      # Documentation (this file)
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

## Quick Start

### Step 1: Clone QE Source

Run the refresh script to clean any existing QE sources and clone the latest QE release:

```powershell
.\scripts\refresh_qe_source.ps1
```

This script will:
- Safely remove QE upstream sources and build artifacts (preserving your code in `scripts/`, `docs/`, etc.)
- Clone the newest QE release from git into `upstream/qe/`
- Initialize all git submodules recursively

**Dry Run Mode:** To see what would be deleted without actually doing it:
```powershell
.\scripts\refresh_qe_source.ps1 -DryRun
```

**Specific Version:** To clone a specific QE version:
```powershell
.\scripts\refresh_qe_source.ps1 -QeVersion qe-7.5
```

### Step 2: Build QE

Run the build script:

```powershell
.\scripts\build_qe_win_oneapi.ps1
```

For a serial build (no MPI):
```powershell
.\scripts\build_qe_win_oneapi.ps1 -NoMpi
```

**Custom QE Source Location:** If QE is cloned to a different location:
```powershell
.\scripts\build_qe_win_oneapi.ps1 -QeSourceDir "path/to/qe"
```

## Build Process

The build script performs these steps:

1. **Verifies QE source and submodules** - Ensures QE is properly cloned with submodules initialized
2. **Discovers Intel oneAPI** - Finds `setvars.bat` to initialize compiler environment
3. **Checks Visual Studio** - Verifies MSVC linker is available
4. **Locates CMake and Ninja** - Finds build tools
5. **Configures CMake** - Sets up build with Intel compilers and MKL
6. **Builds QE** - Compiles `pw.x`, `ph.x`, `pp.x` executables

## CMake/Windows workarounds we had to apply

These changes are already in the repo; keep them in mind for future updates:

- **Preprocessing without shells:** `qe_prepare_fortran_include` (in `cmake/qeHelpers.cmake`) now copies pure Fortran headers and preprocesses the rest via a CMake `-P` helper using `OUTPUT_FILE` (no `>` redirection). It detects `#` lines and picks compiler-specific flags: Intel `-fpp -E -P`, GNU `-cpp -E -P`, NVHPC `-Mpreprocess -E`. Lists are semicolon-escaped when passed with `-D...:STRING=` to avoid CMake argument splitting.
- **LAXlib headers:** `LAXlib/CMakeLists.txt` now uses `qe_prepare_fortran_include` for all `.h -> .fh` generation, eliminating the previous allowlist and shell redirection that broke on Windows.
- **C preprocessor selection:** `QE_CPP` is set explicitly to `icx` in the build script so `qe_preprocess_source` has a working preprocessor on Windows.
- **Git info without pipes:** `cmake/GitInfo.cmake` generates `git-rev.h` via pure `execute_process` + `file(WRITE)`, no `sed`/pipes. Falls back to `UNKNOWN` if git is absent.
- **Submodules must be from git:** QE must be cloned (not zipped) so submodule content and git metadata exist; otherwise CMake will fail. The `refresh_qe_source.ps1` script handles cloning and submodules.
- **Tooling expectations:** Ninja is the generator; VS2022 is required for the MSVC linker; oneAPI `setvars.bat` must be run in the same `cmd.exe` as CMake configure/build (handled by the build script).
- **Internal libs to reduce friction:** We use internal FFTW (`QE_FFTW_VENDOR=Internal`) and internal wannier90/mbd/devxlib to avoid external `find_package` or submodule cloning during configure.

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

- QE upstream source is not committed; Windows fixes live as patches under `patches/`.
- `scripts/refresh_qe_source.ps1` applies these patches automatically (opt-out: `-NoPatch`).
- Manual apply (if you clone QE yourself):
  ```powershell
  .\scripts\apply_qe_patches.ps1 -QeDir upstream/qe -PatchDir patches
  ```
- Patch files:
  - `qe-win-cmake-generation.patch` – CMake/build helpers, Fortran include preprocessing, git-rev generation
  - `qe-win-c-portability.patch` – Windows portability fixes in C sources
  - `qe-win-submodules.patch` – submodule pointer updates (external/d3q, external/devxlib)
- Always clone QE via git with submodules:
  ```powershell
  git clone --recursive https://github.com/QEF/q-e.git
  ```
- To review what’s patched:
  ```powershell
  cd upstream/qe
  git status
  git diff --stat
  ```

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
The script copies QE `.exe` files plus Intel OpenMP, MKL, and Fortran runtime DLLs from your oneAPI install so the binaries run on machines without the full toolchain.

If dumpbin.exe is available (Visual Studio tools), the script will:
- Parse dependencies of each staged `.exe` via `dumpbin /DEPENDENTS`
- Copy the required non-system DLLs from oneAPI redistributable/bin/lib paths (recursively, depth-limited)
- Fall back to pattern-based copies and emit optional-missing notices (`mkl_intel_lp64*.dll`, `libintlc*.dll`)
- Probe the staged binary (default `pw.x.exe -help`) and surface missing-DLL errors if any

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
  - `external/devxlib/src/timer.c` uses `GetSystemTimeAsFileTime`.  
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

## Additional Resources

- QE Official Documentation: https://www.quantum-espresso.org/
- Intel oneAPI Documentation: https://www.intel.com/content/www/us/en/developer/tools/oneapi/documentation.html
- QE GitHub Repository: https://github.com/QEF/q-e


