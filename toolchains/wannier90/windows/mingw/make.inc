# Wannier90 v3.1.0 make.inc for Windows (MinGW-w64 + gfortran)
# Based on make.inc.gfort template
# IMPORTANT: Uses MinGW-w64 gfortran (NOT Intel ifort)
# Serial build by default (MPI disabled)

# Compiler
F90 = gfortran

# Compiler flags
FCOPTS = -O3
LDOPTS = 

# Libraries
# Use OpenBLAS for BLAS/LAPACK on Windows (MinGW-w64)
# OpenBLAS includes both BLAS and LAPACK, so only link -lopenblas
LIBS = -lopenblas

# MPI support (commented out by default)
# To enable parallel postw90:
# COMMS = mpi
# MPIF90 = mpif90

# Installation (optional, used by 'make install')
PREFIX = C:\wannier90
BINDIR = $(PREFIX)\bin

# Utility programs
TAR = tar

