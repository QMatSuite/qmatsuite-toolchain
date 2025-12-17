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
- Replaces hardcoded `.x` suffixes with `${EXE_EXT}` in Makefiles (`Makefile` and `src/Makefile.2`)
- Updates all executable references to use the configurable suffix
- Allows building executables without the `.x` suffix by setting `EXE_EXT=` when building

**What the patch changes:**
- In `Makefile`: Adds `EXE_EXT ?= .x` and replaces hardcoded `.x` in install and clean targets
- In `src/Makefile.2`: Replaces all hardcoded `.x` suffixes in build rules with `${EXE_EXT}`

**Current Status**: This patch is **NOT used** in the CI workflow (`.github/workflows/wannier90-windows-mingw.yml`).

**Why it's not used:**
The patch is unnecessary for our build process because:
1. We build with the default `.x` suffix (which is what Wannier90 builds by default without the patch)
2. The packaging step already handles Windows executable naming by:
   - Finding all `.x` executables: `find "$W90_SRC_DIR" -maxdepth 1 -name "*.x"`
   - Copying them to the bin directory
   - Renaming them to `.x.exe` for Windows: `mv "$f" "${f}.exe"`

Since we're not using `EXE_EXT=` to build without the suffix, and we rename executables in the packaging step anyway, the patch provides no benefit for our workflow.

**Source**: [https://raw.githubusercontent.com/msys2/MINGW-packages/073221b97574014757ce5231e654f085cd6727dd/mingw-w64-wannier90/001-fix-build-on-mingw.patch](https://raw.githubusercontent.com/msys2/MINGW-packages/073221b97574014757ce5231e654f085cd6727dd/mingw-w64-wannier90/001-fix-build-on-mingw.patch)

## Usage

If you want to use this patch (e.g., to build executables without the `.x` suffix), you can apply it manually:

```bash
cd "$W90_SRC_DIR"
patch -Nbp1 -i "$GITHUB_WORKSPACE/toolchains/wannier90/windows/mingw/patches/001-fix-build-on-mingw.patch"
```

Then build with `EXE_EXT=` to create executables without the `.x` suffix:
```bash
make EXE_EXT= wannier90 postw90
```

## Updating Patches

To update patches from the MSYS2 repository, check the latest version at:
https://github.com/msys2/MINGW-packages/tree/master/mingw-w64-wannier90

