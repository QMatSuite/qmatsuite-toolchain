# Wannier90 v3.1.0 make.inc for Linux (GNU Fortran + OpenBLAS)
# Based on make.inc.gfort template

# Compiler
F90 = gfortran

# Compiler flags
F90FLAGS = -O3 -g -fPIC -fmax-errors=1

# Preprocessor
CPP = cpp -P -C -traditional
CPPFLAGS = 

# Linker
LDFLAGS = 

# Libraries
# Use OpenBLAS for BLAS/LAPACK on Linux
LIBS = -lopenblas -llapack

# Parallel build support (optional)
MPIF90 = mpif90
MPILIBS = 

# Installation (optional, used by 'make install')
PREFIX = /usr/local
BINDIR = $(PREFIX)/bin

# Utility programs
TAR = tar

