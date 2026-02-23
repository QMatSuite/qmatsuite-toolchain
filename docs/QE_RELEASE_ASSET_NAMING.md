# QE Release Asset Naming Convention

This document defines the release tag and asset filenames expected by downstream tooling (including main repo installer resolution).

## Windows (from `qe-windows-oneapi-msmpi-release.yml`)

- Tag pattern: `qe-7.5-win-oneapi-msmpi-YYYYMMDD-SHORTSHA`
- Asset: `qe-7.5-win-oneapi-msmpi.zip`
- Alternate asset (libxc build): `qe-7.5-win-oneapi-msmpi-libxc.zip`

## macOS arm64 (from `qe-macos-release.yml`)

- Tag pattern: `qe-7.5-macos-arm64-openmp-YYYYMMDD-SHORTSHA`
- Asset: `qe-7.5-macos-arm64-openmp.zip`
- Optional non-default build-mode tag suffix: `-pw` (only when `build_mode=pw`)

## Linux x64 (future)

- Tag pattern: `qe-7.5-linux-x64-openmp-YYYYMMDD-SHORTSHA`
- Asset: `qe-7.5-linux-x64-openmp.zip`
