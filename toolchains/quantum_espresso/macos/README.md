# QE on macOS

macOS-specific build notes and scripts for Quantum ESPRESSO.

Current CI strategy (see `.github/workflows/qe-linux-macos.yml`):

- Runner: `macos-14`
- Dependencies installed via Homebrew
- BLAS/LAPACK from Apple's **Accelerate** framework
- Build at least `pw.x` using `make pw`

