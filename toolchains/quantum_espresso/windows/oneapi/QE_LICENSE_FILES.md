# Quantum ESPRESSO License Files for Redistribution

This document lists all license files found in the QE source tree that should be included in the redistributed binary distribution.

## Required License Files (Main QE)

1. **`License`** (root)
   - Type: GNU General Public License v2 (GPL-2.0)
   - Location: `upstream/qe/License`
   - Status: **REQUIRED** - This is the main QE license

## Component-Specific License Files (Should Include)

2. **`EPW/License`**
   - Type: GNU General Public License v2 (GPL-2.0)
   - Location: `upstream/qe/EPW/License`
   - Status: **SHOULD INCLUDE** - EPW module has its own license file

3. **`external/wannier90/LICENSE`**
   - Type: GNU General Public License v2 (GPL-2.0)
   - Location: `upstream/qe/external/wannier90/LICENSE`
   - Status: **SHOULD INCLUDE** - Wannier90 is a bundled external library

4. **`external/lapack/LICENSE`**
   - Type: BSD-style (University of Tennessee/California/Colorado)
   - Location: `upstream/qe/external/lapack/LICENSE`
   - Status: **SHOULD INCLUDE** - LAPACK is bundled and has different license terms

5. **`external/fox/LICENSE`**
   - Type: BSD-style (FoX XML library)
   - Location: `upstream/qe/external/fox/LICENSE`
   - Status: **SHOULD INCLUDE** - FoX is bundled and has different license terms

6. **`external/mbd/LICENSE`**
   - Type: Mozilla Public License v2.0 (MPL-2.0)
   - Location: `upstream/qe/external/mbd/LICENSE`
   - Status: **SHOULD INCLUDE** - MBD has different license (MPL-2.0)

7. **`external/qe-gipaw/LICENSE`**
   - Type: GNU General Public License v2 (GPL-2.0)
   - Location: `upstream/qe/external/qe-gipaw/LICENSE`
   - Status: **SHOULD INCLUDE** - qe-gipaw is a bundled external library

8. **`external/d3q/License`**
   - Type: GPL v2 (with CeCILL note)
   - Location: `upstream/qe/external/d3q/License`
   - Status: **SHOULD INCLUDE** - d3q is bundled and has licensing notes

## Optional Documentation Files

9. **`README.md`** (root)
   - Location: `upstream/qe/README.md`
   - Status: **OPTIONAL** - Main README, mentions license badge but not required for legal compliance

10. **`CONTRIBUTING.md`** (root)
    - Location: `upstream/qe/CONTRIBUTING.md`
    - Status: **OPTIONAL** - Contributing guidelines, not required for legal compliance

## GUI Component Licenses (Not Required for Binary Distribution)

- `GUI/PWgui/COPYING` - GUI component, not part of binary distribution
- `GUI/Guib/COPYING` - GUI component, not part of binary distribution
- Various `license.terms` files in GUI directories - GUI components only

## Recommendation

**For binary redistribution, include:**
1. Main `License` file (REQUIRED)
2. All external library licenses (items 2-8 above) for complete legal compliance

**File structure in dist/licenses/quantum-espresso/:**
```
dist/licenses/quantum-espresso/
├── License                    (main QE license)
├── EPW-License                (EPW module)
├── wannier90-LICENSE          (external/wannier90)
├── lapack-LICENSE             (external/lapack)
├── fox-LICENSE                (external/fox)
├── mbd-LICENSE                (external/mbd)
├── qe-gipaw-LICENSE           (external/qe-gipaw)
└── d3q-License                (external/d3q)
```

Note: Files are renamed to avoid conflicts (e.g., multiple "LICENSE" files would overwrite each other).

