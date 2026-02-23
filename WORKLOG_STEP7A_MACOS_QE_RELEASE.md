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
