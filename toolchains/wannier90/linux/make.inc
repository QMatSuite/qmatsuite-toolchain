# Wannier90 v3.1.0 make.inc for Linux (GNU Fortran)
# Based on make.inc.gfort template
# Serial build by default (MPI disabled)

# Compiler
F90 = gfortran

# Compiler flags
FCOPTS = -O3
LDOPTS = 

# Libraries
# Use system BLAS/LAPACK
LIBS = -llapack -lblas

# MPI support (commented out by default)
# To enable parallel postw90:
# COMMS = mpi
# MPIF90 = mpif90

# Installation (optional, used by 'make install')
PREFIX = /usr/local
BINDIR = $(PREFIX)/bin

# Utility programs
TAR = tar

