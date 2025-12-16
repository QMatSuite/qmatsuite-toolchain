# Wannier90 v3.1.0 make.inc for Windows (Intel Fortran Classic)
# Based on make.inc.ifort template
# IMPORTANT: Uses Intel Fortran Classic (ifort), NOT MinGW gfortran

# Compiler
F90 = ifort

# Compiler flags
# Use -fpp for Intel Fortran preprocessor (equivalent to -cpp for gfortran)
F90FLAGS = -O3 -fpp

# Preprocessor
CPP = ifort -fpp -E
CPPFLAGS = 

# Linker
LDFLAGS = 

# Libraries
# Use Intel MKL if available, otherwise use system BLAS/LAPACK
# MKL is typically available with Intel oneAPI
LIBS = -mkl

# Parallel build support (optional)
MPIF90 = mpiifort
MPILIBS = 

# Installation (optional, used by 'make install')
PREFIX = C:\wannier90
BINDIR = $(PREFIX)\bin

# Utility programs
TAR = tar

