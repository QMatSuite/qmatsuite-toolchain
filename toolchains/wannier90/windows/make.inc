# Wannier90 v3.1.0 make.inc for Windows (Intel Fortran Classic)
# Based on make.inc.ifort template
# IMPORTANT: Uses Intel Fortran Classic (ifort), NOT MinGW gfortran
# Serial build by default (MPI disabled)

# Compiler
F90 = ifort

# Compiler flags
# -fpp: use Fortran preprocessor (for consistency with QE Windows builds)
FCOPTS = -O2 -fpp
LDOPTS = -O2

# Libraries
# Use Intel MKL (sequential mode to avoid OpenMP conflicts)
# If -mkl=sequential fails, fall back to explicit MKL libs using MKLROOT
LIBS = -mkl=sequential

# MPI support (commented out by default)
# To enable parallel postw90:
# COMMS = mpi
# MPIF90 = mpiifort

# Installation (optional, used by 'make install')
PREFIX = C:\wannier90
BINDIR = $(PREFIX)\bin

# Utility programs
TAR = tar

