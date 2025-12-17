# Quantum ESPRESSO — Windows oneAPI (ifx/icx + CMake/Ninja)

This folder packages the Windows native oneAPI workflow for Quantum ESPRESSO (QE). It clones QE, applies Windows-focused patches, builds with ifx/icx via CMake/Ninja, and stages runnable binaries.

## One-click flow (run from repo root)

```powershell
.\toolchains\quantum_espresso\windows\oneapi\scripts\refresh_qe_source.ps1
.\toolchains\quantum_espresso\windows\oneapi\scripts\build_qe_win_oneapi.ps1 -NoMpi
.\toolchains\quantum_espresso\windows\oneapi\scripts\stage_qe_windows.ps1
```

- QE is cloned to `toolchains/quantum_espresso/windows/oneapi/upstream/qe`
- Build output lives in `.../upstream/qe/build-win-oneapi`
- Staged binaries land in `.../dist/win-oneapi`

## Patches (applied automatically by refresh script)

- `qe-win-cmake-generation.patch` — CMake helpers, Fortran include preprocessing without shell redirection, git-rev header generation
- `qe-win-c-portability.patch` — Windows portability fixes in QE C sources
- `qe-win-submodules.patch` — submodule pointer updates (e.g., external/d3q, external/devxlib)

Manual apply if needed:
```powershell
.\toolchains\quantum_espresso\windows\oneapi\scripts\apply_qe_patches.ps1 `
  -QeDir toolchains/quantum_espresso/windows/oneapi/upstream/qe `
  -PatchDir toolchains/quantum_espresso/windows/oneapi/patches
```

## Why git clone (not zip)

QE must be obtained via `git clone --recursive` so submodules and git metadata exist; zip/tar archives break submodules and git-derived version info. The refresh script handles clone + submodules and then applies the local patches.

## Build notes

- Toolchain: Intel oneAPI (ifx/icx) + VS2022 linker, CMake, Ninja.
- MPI: pass `-NoMpi` for serial; omit to enable MPI if available.
- CMake flags used (key ones): `-G Ninja`, `-DQE_CPP=<icx>`, `-DQE_LAPACK_INTERNAL=OFF`, `-DQE_WANNIER90_INTERNAL=ON`, `-DQE_MBD_INTERNAL=ON`, `-DQE_DEVICEXLIB_INTERNAL=ON`, `-DQE_FFTW_VENDOR=Internal`.

## Staging distributables

`scripts/stage_qe_windows.ps1` copies built `.exe` and Intel runtime DLLs into `dist/win-oneapi/`. If `dumpbin.exe` is available it scans dependencies; otherwise it falls back to pattern-based copying. Optional flags: `-Only`, `-Clean`, `-VerboseDeps`.

## Layout

```
toolchains/quantum_espresso/windows/oneapi/
  README.md
  scripts/
    refresh_qe_source.ps1
    apply_qe_patches.ps1
    build_qe_win_oneapi.ps1
    stage_qe_windows.ps1
  patches/
    qe-win-cmake-generation.patch
    qe-win-c-portability.patch
    qe-win-submodules.patch
  tests/   (placeholder)
```

