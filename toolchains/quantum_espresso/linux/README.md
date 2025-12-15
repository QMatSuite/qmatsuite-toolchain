# QE on Linux

Linux-specific build notes and scripts for Quantum ESPRESSO.

Current CI strategy (see `.github/workflows/qe-linux-macos.yml`):

- Runner: `ubuntu-latest`
- Dependencies installed via `apt`
- QE 7.5 source downloaded from the official GitLab archive
- Configure with system BLAS/LAPACK (OpenBLAS)
- Build at least `pw.x` using `make pw`

