# Fix for MinGW gfortran i2c() Empty String Bug

## Root Cause

On MinGW gfortran (MSYS2 UCRT64), the `i2c()` function in `xmltools.f90` returns an empty string when converting integers to character strings. This is due to a compiler-specific bug with allocatable character functions.

The `i2c()` function was declared as:
```fortran
function i2c (i) result (c)
  character(len=:), allocatable :: c
  ...
end function i2c
```

When used in tag construction like `'PP_BETA.'//i2c(nb)`, the empty result causes tags to become `PP_BETA.` instead of `PP_BETA.1`, `PP_BETA.2`, etc. This leads to:
1. Tag mismatch (UPF file has `PP_BETA.1`, code searches for `PP_BETA.`)
2. Empty arrays returned from `xmlr_readtag()`
3. Zero-length arrays passed to `simpson()` integration
4. Runtime crash: "Index '1' of array 'func' above upper bound of 0"

## Fix Approach

### 1. Fixed `i2c()` to return fixed-length CHARACTER
- Changed from `character(len=:), allocatable :: c` to `character(len=32) :: c`
- Added safety check to detect and abort if conversion results in empty string
- Uses `I0` format instead of `i11` for more efficient formatting

### 2. Replaced `i2c()` usage in tag construction
- All tag constructions now use direct `WRITE(nbuf,'(I0)')` to fixed buffer
- Avoids dependency on `i2c()` entirely for critical tag construction
- Added safety checks to ensure tags are non-empty after construction

### 3. Added optional diagnostics
- Debug output controlled by `QE_UPF_DEBUG` environment variable
- When set, prints constructed tag strings for `PP_BETA.x` and `PP_CHI.x`
- Helps verify correct tag construction during debugging

### 4. Added robust array size checks
- After reading each `PP_BETA.x` and `PP_CHI.x` block, verify array size > 0
- If empty, produce controlled error with context (which tag, which index)
- Prevents silent failures and downstream crashes in `simpson()`

## Files Modified

- `upflib/xmltools.f90`: Fixed `i2c()` function
- `upflib/read_upf_new.f90`: Replaced all `i2c()` usages in tag construction

## Testing

1. **MinGW Build**: Rebuild QE with this patch applied
2. **Run Test**: Use `scf-cg.in` with `Si.pz-vbc.UPF`
3. **Verify**: 
   - No crash in `simpson()`
   - Correct `angular_momentum` values parsed (0 and 1, not both 0)
   - Arrays are non-empty before reaching `simpson()`

## Debugging

To enable diagnostic output:
```bash
export QE_UPF_DEBUG=1
pw.x -in scf-cg.in > scf-cg.out 2>&1
```

Look for lines like:
```
DEBUG read_pp_nonlocal: nb=1, constructed tag=[PP_BETA.1]
DEBUG read_pp_nonlocal: nb=2, constructed tag=[PP_BETA.2]
DEBUG read_pp_pswfc: nw=1, constructed tag=[PP_CHI.1]
DEBUG read_pp_pswfc: nw=2, constructed tag=[PP_CHI.2]
```

## Compatibility

- **Intel oneAPI**: No changes to behavior (already worked correctly)
- **MinGW gfortran**: Fixes the empty string bug
- **Other compilers**: Should work correctly with fixed-length CHARACTER

## Why MinGW Triggered It But Intel Didn't

MinGW gfortran has a known issue with allocatable character functions in certain contexts (likely related to internal file I/O and string concatenation). Intel oneAPI's Fortran compiler handles allocatable characters more robustly in this scenario. The fix (using fixed-length CHARACTER) works correctly on both compilers.

