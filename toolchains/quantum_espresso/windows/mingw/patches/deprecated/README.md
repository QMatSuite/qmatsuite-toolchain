# Deprecated Patches

This directory contains intermediate patches that were produced during debugging and have been replaced by a single unified patch.

## Replaced Patches

- **005-fix-xmltools-eof-rewind.patch** - Defensive EOF rewind fix
- **006-fix-empty-tag-robust-mingw.patch** - `capitalize_if_v2` fix
- **006-fix-empty-tag-in-read-pp-header-v2.patch** - Alternative empty tag fix (superseded)
- **007-fix-i2c-mingw-empty-string.patch** - `i2c()` fix and tag construction updates
- **007-fix-i2c-mingw-empty-string-README.md** - Documentation for patch 007

## Replacement

All fixes from patches 005-007 have been consolidated into a single authoritative patch:

**`005-mingw-upf-xml-fixes.patch`**

This unified patch contains all final, working fixes related to:
- UPF v2 XML parsing issues (empty tag, attribute parsing, MinGW gfortran quirks)
- `upflib/xmltools.f90` fixes (i2c, empty tag guard, defensive EOF rewind)
- `upflib/read_upf_new.f90` fixes (capitalize_if_v2_into, i2c replacement, robust checks)
- `upflib/write_upf_new.f90` fixes (capitalize_if_v2_into for consistency)

## Why Deprecated

These patches were created incrementally during the debugging process. The unified patch:
- Consolidates all fixes into a single authoritative patch
- Ensures consistency across all three files
- Provides a single, authoritative source of truth
- Matches the known-good working state exactly (including debug instrumentation)

## Usage

**Do not apply these deprecated patches.** Use `005-mingw-upf-xml-fixes.patch` instead.

