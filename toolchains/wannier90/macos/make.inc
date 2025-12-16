# Wannier90 v3.1.0 make.inc for macOS (GNU Fortran + Apple Accelerate)
# Based on make.inc.macosx.macports template
# IMPORTANT: Uses Apple Accelerate framework with -ff2c flag (NOT OpenBLAS)
# Serial build by default (MPI disabled due to -ff2c/OpenMPI wrapper incompatibilities)

# Compiler
F90 = gfortran
CC = cc

# Compiler flags
# -ff2c is required for Accelerate framework compatibility
FCOPTS = -O2 -ff2c
LDOPTS = -O2 -ff2c

# Libraries
# MUST use Apple Accelerate framework for BLAS/LAPACK on macOS
# Do NOT use OpenBLAS as the default
LIBS = -framework Accelerate

# MPI support (commented out by default)
# WARNING: -ff2c can conflict with OpenMPI wrappers
# To enable MPI, remove -ff2c from FCOPTS and LDOPTS first
# COMMS = mpi
# MPIF90 = mpif90

# Installation (optional, used by 'make install')
PREFIX = /usr/local
BINDIR = $(PREFIX)/bin

# Utility programs
TAR = tar

