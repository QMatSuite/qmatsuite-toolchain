# QE Windows Patch Notes

This repo carries small Windows-focused patches for the upstream QE sources.  
Apply them after cloning QE (e.g., after running `scripts/refresh_qe_source.ps1`) using `scripts/apply_qe_patches.ps1`.

## Patch list

- `qe-win-cmake-generation.patch`  
  Makes QE CMake generation Windows-friendly (avoids shell pipelines and uses CMake scripting for git info and Fortran includes).

- `qe-win-c-portability.patch`  
  Adds Windows stubs/compatibility for socket and math code.

- `qe-win-devxlib-timer.patch`  
  Replaces the POSIX-only `sys/times.h` usage in `external/devxlib/src/timer.c` with a Windows `QueryPerformanceCounter` implementation so DeviceXlib builds on Windows.

## How to apply

From the repo root:

```powershell
.\scripts\apply_qe_patches.ps1
```

The script is idempotent: it skips patches already applied and stops if a patch cannot be applied cleanly.

## Workflow reminder

- `scripts\refresh_qe_source.ps1` only cleans and clones QE (no patching).  
- After cloning, run `scripts\apply_qe_patches.ps1` to apply the Windows patches above.  
- Then build with `scripts\build_qe_win_oneapi.ps1` (add `-NoMpi` if you want serial builds).

## Recent fix: devxlib timer patch

- The `qe-win-devxlib-timer.patch` had path/EOL issues; it has been rebuilt and now applies cleanly against `external/devxlib/src/timer.c`.  
- If you previously saw “corrupt patch” errors, sync the updated patch file and re-run `scripts\apply_qe_patches.ps1`.

