# QE on Windows

Windows-specific build notes and scripts for Quantum ESPRESSO.

Current CI strategy (see `.github/workflows/qe-windows-intel.yml`):

- Runner: `windows-latest`
- Install **Intel oneAPI HPC toolkit** (Fortran compiler + MKL) via
  `modflowpy/install-intelfortran-action`.
- Install **MS-MPI runtime + SDK** using Chocolatey (`msmpi`, `msmpisdk`).
- Use `ifort` as the Fortran compiler, with MKL for BLAS/LAPACK.
- Build at least `pw.x` (OpenMP build, MPI is present for future use).

The end goal for the distributed Windows package is a layout similar to classic
precompiled QE distributions:

- Many `*.exe` programs (`pw.exe`, `ph.exe`, `bands.exe`, etc.)
- MKL runtime DLLs (`mkl_*.dll`)
- MPI executables (`mpiexec.exe`, `smpd.exe`) and `msmpi.dll`

