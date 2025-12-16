# Wannier90 Windows MinGW Patches

This directory contains patches for building Wannier90 on Windows using MinGW-w64.

## Source

**Patches in this directory are sourced from the MSYS2 MINGW-packages repository:**

**Source Repository**: [https://github.com/msys2/MINGW-packages/tree/master/mingw-w64-wannier90](https://github.com/msys2/MINGW-packages/tree/master/mingw-w64-wannier90)

These patches are maintained by the MSYS2 project and are used in their official Wannier90 package for MinGW-w64.

## Patches

### 001-fix-build-on-mingw.patch

This patch fixes executable suffix handling in Wannier90 Makefiles for Windows/MinGW. It:

- Adds an `EXE_EXT` variable (defaults to `.x`) that can be set to empty on Windows
- Replaces hardcoded `.x` suffixes with `${EXE_EXT}` in Makefiles
- Allows building executables without the `.x` suffix on Windows

**Current Status**: This patch is **currently applied** in the CI workflow (`.github/workflows/wannier90-windows-mingw.yml`) before building. The patch makes the executable suffix configurable via the `EXE_EXT` variable, which defaults to `.x`. The current build process uses the default `.x` suffix, so executables are built as `wannier90.x`, `postw90.x`, etc.

**Source**: [https://raw.githubusercontent.com/msys2/MINGW-packages/073221b97574014757ce5231e654f085cd6727dd/mingw-w64-wannier90/001-fix-build-on-mingw.patch](https://raw.githubusercontent.com/msys2/MINGW-packages/073221b97574014757ce5231e654f085cd6727dd/mingw-w64-wannier90/001-fix-build-on-mingw.patch)

## Usage

This patch is automatically applied by the CI workflow (`.github/workflows/wannier90-windows-mingw.yml`) before building.

The patch is applied with:
```bash
cd "$W90_SRC_DIR"
patch -Np1 -i "$GITHUB_WORKSPACE/toolchains/wannier90/windows/patches/001-fix-build-on-mingw.patch"
```

By default, the build uses `EXE_EXT=.x` (the default value), so executables are built with the `.x` suffix. To build without the suffix, you can set `EXE_EXT=` when building:
```bash
make EXE_EXT= wannier90 postw90
```

## Updating Patches

To update patches from the MSYS2 repository, check the latest version at:
https://github.com/msys2/MINGW-packages/tree/master/mingw-w64-wannier90

