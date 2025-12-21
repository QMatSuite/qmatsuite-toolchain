# Third-Party License Files for Redistribution

This document lists all third-party license files that should be included in the binary distribution (besides Quantum ESPRESSO licenses).

## Intel oneAPI Runtime Libraries

### Found License Files:

1. **Main oneAPI License**
   - Location: `C:\Program Files (x86)\Intel\oneAPI\licensing\2025.3\license.htm`
   - Type: HTML license file
   - Status: **SHOULD INCLUDE** - Main oneAPI license

2. **Installer License**
   - Location: `C:\Program Files (x86)\Intel\oneAPI\Installer\license.txt`
   - Type: Text license file
   - Status: **OPTIONAL** - Installer license, may not be needed for runtime

3. **Component-Specific Licenses** (found in various components):
   - `C:\Program Files (x86)\Intel\oneAPI\advisor\2025.4\licensing\license.txt`
   - `C:\Program Files (x86)\Intel\oneAPI\vtune\2025.7\licensing\license.txt`
   - `C:\Program Files (x86)\Intel\oneAPI\2025.3\include\stb\LICENSE.txt` (third-party library)

### Recommended Search Strategy:

Search for license files in:
- `C:\Program Files (x86)\Intel\oneAPI\licensing\` - Main license directory
- `C:\Program Files (x86)\Intel\oneAPI\compiler\*\licensing\` - Compiler licenses
- `C:\Program Files (x86)\Intel\oneAPI\mkl\*\licensing\` - MKL licenses
- `C:\Program Files (x86)\Intel\oneAPI\*\redist\*\` - Redistributable licenses

**Patterns to search:**
- `LICENSE*` (case-insensitive)
- `THIRD-PARTY*`
- `NOTICE*`
- `license.txt`
- `license.htm`

## Microsoft MS-MPI

### Found License Files:

1. **MS-MPI SDK License**
   - Location: `C:\Program Files (x86)\Microsoft SDKs\MPI\License`
   - Type: License file (no extension)
   - Status: **SHOULD INCLUDE** - MS-MPI SDK license

2. **MS-MPI Runtime License**
   - Location: `C:\Program Files\Microsoft MPI\License`
   - Type: License file (no extension)
   - Status: **SHOULD INCLUDE** - MS-MPI runtime license

### Additional Search Locations:
- `C:\Program Files (x86)\Microsoft SDKs\MPI\License.txt`
- `C:\Program Files (x86)\Microsoft SDKs\MPI\LICENSE.txt`
- `C:\Program Files (x86)\Microsoft SDKs\MPI\EULA.txt`
- `C:\Program Files\Microsoft MPI\License.txt`

## Microsoft Visual C++ Runtime

### Typical License Locations:

1. **VC++ Redistributable License**
   - Typical location: `C:\Program Files (x86)\Microsoft Visual Studio\*\VC\Redist\*\license.txt`
   - Type: Text license file
   - Status: **SHOULD INCLUDE** - VC++ runtime redistributable license

2. **Alternative Locations:**
   - `C:\Program Files\Microsoft Visual Studio\*\VC\Redist\*\license.txt`
   - `C:\Program Files (x86)\Microsoft Visual C++ Redistributable\license.txt`

### Search Strategy:

Search for `license.txt` in:
- Visual Studio installation directories under `VC\Redist\`
- Microsoft Visual C++ Redistributable directories

## Recommended Implementation

### Directory Structure:
```
dist/licenses/
├── quantum-espresso/
│   ├── License                    (main QE license)
│   ├── EPW/
│   │   └── License
│   └── external/
│       ├── wannier90/
│       │   └── LICENSE
│       ├── lapack/
│       │   └── LICENSE
│       ├── fox/
│       │   └── LICENSE
│       ├── mbd/
│       │   └── LICENSE
│       ├── qe-gipaw/
│       │   └── LICENSE
│       └── d3q/
│           └── License
├── intel-oneapi/
│   ├── license.htm                (main oneAPI license)
│   ├── license.txt                (if found)
│   └── [component-specific licenses]
├── microsoft-msmpi/
│   └── License                    (MS-MPI license)
└── microsoft-vs-runtime/
    └── license.txt                (VC++ runtime license)
```

## Files to Include (Summary)

### Intel oneAPI:
- ✅ Main license from `licensing\` directory
- ✅ Any redistributable-specific licenses
- ⚠️ Component licenses (compiler, MKL) - include if found in redist directories

### Microsoft MS-MPI:
- ✅ `C:\Program Files (x86)\Microsoft SDKs\MPI\License`
- ✅ `C:\Program Files\Microsoft MPI\License` (if different from SDK)

### Microsoft VC++ Runtime:
- ✅ `license.txt` from VC++ Redistributable directories

## Notes

- Intel oneAPI licenses may be in HTML format (`license.htm`) - include as-is
- MS-MPI license files may not have extensions - preserve original names
- VC++ runtime license is typically a single `license.txt` file
- Some licenses may be version-specific (e.g., `2025.3\license.htm`) - include the most recent or all found

