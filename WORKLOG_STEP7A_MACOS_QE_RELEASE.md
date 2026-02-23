# WORKLOG STEP7A: macOS QE Release Workflow Upgrade

## 2026-02-23

### 1) Context Read-in (Mandatory)
- Read `.github/workflows/qe-linux-macos.yml`.
- Read `.github/workflows/qe-windows-oneapi-msmpi-release.yml` as template/gold standard.
- Inspected QE toolchain directory structure under `toolchains/quantum_espresso`.
- Verified platform directories:
  - `toolchains/quantum_espresso/macos/README.md`
  - `toolchains/quantum_espresso/linux/README.md`
- Read Windows staging script:
  - `toolchains/quantum_espresso/windows/oneapi/scripts/stage_qe_windows.ps1`

### 2) Findings
- Existing Linux/macOS workflow is CI-oriented (push/PR), builds `pw.x` only, uploads artifact only.
- Windows release workflow includes release-grade flow: manual dispatch, packaging, checksums, attestation, release tag handling, GitHub Release upload.
- Existing repository has no dedicated macOS QE release workflow file.

### 3) Plan
- Add new workflow: `.github/workflows/qe-macos-release.yml`.
- Keep existing `.github/workflows/qe-linux-macos.yml` intact for CI.
- Implement manual release build with:
  - `build_mode` (`pw` / `all`)
  - release tag auto-generation if empty
  - optional draft release
  - full staging (`.x` executables) to `bin/`
  - zip + SHA256 checksums
  - artifact attestation
  - `softprops/action-gh-release@v2` release publishing
- Add documentation for release asset naming for main repo integration.

### 4) Implementation Completed

#### 4.1 New workflow added
- Created `.github/workflows/qe-macos-release.yml`.
- Trigger mode: `workflow_dispatch` only.
- Inputs implemented:
  - `build_mode` (`pw` / `all`)
  - `release_tag` (optional override)
  - `release_draft` (boolean)
  - `sign` (boolean; default `true`) for macOS codesign/notarization

#### 4.2 Build pipeline behavior
- Runner: `macos-14` (Apple Silicon).
- Downloads QE 7.5 source tarball.
- Configures QE for serial build with Accelerate:
  - `CC`, `F90`, `MPIF90`, `CPP`
  - `BLAS_LIBS=-framework Accelerate`
  - `LAPACK_LIBS=-framework Accelerate`
- Build mode:
  - `all` -> `make all -j$(sysctl -n hw.ncpu)`
  - `pw` -> `make pw -j$(sysctl -n hw.ncpu)`
- Validation:
  - always requires `pw.x`
  - for `all`, requires: `pw.x ph.x dos.x bands.x pp.x`

#### 4.3 Staging and packaging
- Stages executables into `$QE_STAGE_DIR/bin`:
  - copies from `QE_SRC_DIR/bin/*.x` first
  - then copies additional executable `*.x` discovered in source tree if not already present
- Creates release zip:
  - `qe-7.5-macos-arm64-openmp.zip`
- Generates checksum file:
  - `checksums.txt` with SHA256

#### 4.4 Release, attestation, and signing alignment
- Adds build provenance attestation (`actions/attest-build-provenance@v2`).
- Creates release via `softprops/action-gh-release@v2` with zip + checksum assets.
- Auto tag generation implemented:
  - default: `qe-7.5-macos-arm64-openmp-YYYYMMDD-SHORTSHA`
  - `build_mode=pw`: `qe-7.5-macos-arm64-openmp-pw-YYYYMMDD-SHORTSHA`
- Aligned to `/Users/hh7465/QMatSuite/.github/workflows/release-macos.yml` for code signing secrets and flow:
  - `APPLE_CERTIFICATE_P12`
  - `APPLE_CERTIFICATE_PASSWORD`
  - `DEVELOPER_NAME`
  - `APPLE_TEAM_ID`
  - `APPLE_NOTARY_ISSUER_ID`
  - `APPLE_NOTARY_KEY_ID`
  - `APPLE_NOTARY_KEY`
- Signing flow:
  - import Apple certificate into temporary keychain
  - `codesign` all staged `*.x` executables
  - notarize the generated zip with `xcrun notarytool`

### 5) Asset Naming Documentation
- Added `docs/QE_RELEASE_ASSET_NAMING.md` with explicit naming patterns for:
  - Windows (current)
  - macOS arm64 (this step)
  - Linux x64 (future)
- Added link entry in `docs/README.md`.

### 6) Verification Status
- Local verification performed:
  - file existence and content checks for new workflow and docs
  - YAML syntax sanity check passed for `.github/workflows/qe-macos-release.yml`
- Remote verification (GitHub Actions) not executable from this environment:
  - `gh auth status` confirms authenticated GitHub CLI session with `workflow` scope.
  - Attempted dispatch:
    - `gh workflow run qe-macos-release.yml -R QMatSuite/qmatsuite-toolchain -f build_mode=all -f release_draft=true -f sign=false`
    - result: `HTTP 404 workflow ... not found on the default branch` (expected before pushing this new workflow file).
  - cannot observe draft release artifacts directly until workflow file is pushed and run.

### 7) Manual Verification Checklist (for repo maintainer)
1. Trigger `.github/workflows/qe-macos-release.yml` with:
   - `build_mode=all`
   - `release_draft=true`
   - `sign=true` (or `false` if secrets are not configured)
2. Confirm job succeeds on `macos-14`.
3. Confirm draft release is created with assets:
   - `qe-7.5-macos-arm64-openmp.zip`
   - `checksums.txt`
4. Download zip and verify required executables present in `bin/`:
   - `pw.x`, `ph.x`, `dos.x`, `bands.x`, `pp.x`
5. Verify checksum:
   - `shasum -a 256 qe-7.5-macos-arm64-openmp.zip`
   - compare with `checksums.txt`
6. Smoke check binaries:
   - `./pw.x --version`
   - `./ph.x --help`
   - `./dos.x --help`

## Step 7A-fix: Portable macOS QE + make all fixes (2026-02-23)

### A) Files updated
- `.github/workflows/qe-macos-release.yml`
- `scripts/bundle_macos_dylibs.sh` (new)

### B) Applied build fixes from main repo tests workflow
- Added PIOUD `-pg` patch step before build.
- Changed build order for `build_mode=all` to:
  - `make pw -j$(sysctl -n hw.ncpu) || make pw -j1`
  - `make all -j1`
- Added wannier90 re-staging from `external/wannier90/*.x` and Mach-O validation checks for critical executables.
- Added FFTW-aware macOS configure flags and `make.inc` FFTW3 define/include patching.

### C) Portability bundling implementation
- Added `scripts/bundle_macos_dylibs.sh`:
  - Collects non-system dylibs recursively (transitive closure)
  - Copies dylibs into `lib/`
  - Rewrites executable load paths to `@executable_path/../lib/...`
  - Rewrites dylib IDs and inter-dylib paths to `@loader_path/...`
  - Re-signs modified binaries ad-hoc for arm64 validity
  - Verifies no `/opt/homebrew` or `/usr/local` references remain

### D) Local verification run (executed)
Local run directory: `/tmp/qe-step7a-fix`

Build status evidence from local log:
- `[INFO] Build pw then all` at line 133
- `[INFO] Stage executables` at line 5986
- `[INFO] otool before bundling` at line 5987
- `[INFO] Bundle dylibs` at line 5994

Interpretation: `make pw` + `make all -j1` completed and reached staging/bundling.

### E) `otool -L` before bundling (`pw.x`)
```text
/tmp/qe-step7a-fix/qe-stage/bin/pw.x:
	/opt/homebrew/opt/fftw/lib/libfftw3.3.dylib (compatibility version 10.0.0, current version 10.10.0)
	/opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.5.dylib (compatibility version 6.0.0, current version 6.0.0)
	/System/Library/Frameworks/Accelerate.framework/Versions/A/Accelerate (compatibility version 1.0.0, current version 4.0.0)
	/opt/homebrew/opt/gcc/lib/gcc/current/libquadmath.0.dylib (compatibility version 1.0.0, current version 1.0.0)
	/usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1351.0.0)
```

### F) `otool -L` after bundling (`pw.x`)
```text
/tmp/qe-step7a-fix/qe-stage/bin/pw.x:
	@executable_path/../lib/libfftw3.3.dylib (compatibility version 10.0.0, current version 10.10.0)
	@executable_path/../lib/libgfortran.5.dylib (compatibility version 6.0.0, current version 6.0.0)
	/System/Library/Frameworks/Accelerate.framework/Versions/A/Accelerate (compatibility version 1.0.0, current version 4.0.0)
	@executable_path/../lib/libquadmath.0.dylib (compatibility version 1.0.0, current version 1.0.0)
	/usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1351.0.0)
```

### G) Clean-PATH portability probe output
Command context: clean environment `PATH=/usr/bin:/bin:/usr/sbin:/sbin`

```text
Program PWSCF v.7.5 starts on 22Feb2026 at 22:59:34

This program is part of the open-source Quantum ESPRESSO suite
for quantum simulation of materials; please cite
```

Result: PASS (no `Library not loaded` errors).

### H) Bundle script summary output
```text
=== Non-system dylibs to bundle ===
/opt/homebrew/opt/fftw/lib/libfftw3.3.dylib
/opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.5.dylib
/opt/homebrew/opt/gcc/lib/gcc/current/libquadmath.0.dylib
=== Total: 3 ===
Bundled: libfftw3.3.dylib
Bundled: libgfortran.5.dylib
Bundled: libquadmath.0.dylib

=== Portability verification ===
All binaries and dylibs are self-contained ✓

=== Bundle summary ===
Executables: 89
Bundled dylibs: 3
Total bin/ size: 393M
Total lib/ size: 3.1M
```

### I) Bundle size summary (recorded)
```text
bin_size=393M
lib_size=3.1M
bin_count=89
lib_count=3
```

### J) Post-fix script correction
- Initial bundling script revision had a transitive-deps loop hazard (reading and appending same file in-loop).
- Fixed by iterating over a snapshot file (`mktemp`) each pass before appending newly discovered deps.

### K) Asset naming confirmation
- Zip asset name remains unchanged: `qe-7.5-macos-arm64-openmp.zip`
- Internal layout now includes both `bin/` and `lib/`.

## Step 7A-cache: Split build/sign and add stage cache (2026-02-23)

### Goal
Avoid recompiling QE when signing/notarization fails and the workflow is rerun, while still forcing recompile when build logic/flags change.

### Implementation
Updated `.github/workflows/qe-macos-release.yml` to two-job flow:

1. `build_stage`
- Builds + stages + bundles QE (or skips build on cache hit).
- Caches staged distribution via `actions/cache@v4`:
  - Path: `.cache/qe-macos-release/stage`
  - Key:
    - `qe-macos-stage-${{ runner.os }}-${{ inputs.build_mode }}-${{ hashFiles('.github/workflows/qe-macos-release.yml', 'scripts/bundle_macos_dylibs.sh') }}`
- Uploads staged artifact:
  - `qe-macos-stage-${{ github.sha }}-${{ inputs.build_mode }}`

2. `sign_package_attest_release`
- Depends on `build_stage` (`needs: build_stage`).
- Downloads staged artifact.
- Runs signing/notarization/package/release only.

### Rerun behavior
- If signing/notarization fails, `Re-run failed jobs` reruns only `sign_package_attest_release` (no compile).
- If full workflow is rerun with same build logic and same `build_mode`, stage cache can short-circuit compile in `build_stage`.
- If build flags/workflow logic or bundling script changes, `hashFiles(...)` changes and cache key invalidates, forcing rebuild.

### Validation
- YAML syntax check passed after change.
