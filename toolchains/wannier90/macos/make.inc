# Wannier90 v3.1.0 make.inc for macOS (GNU Fortran + Apple Accelerate)
# Based on make.inc.macosx.homebrew template
# IMPORTANT: Uses Apple Accelerate framework (NOT OpenBLAS)

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
# MUST use Apple Accelerate framework for BLAS/LAPACK on macOS
# Do NOT use OpenBLAS as the default
LIBS = -framework Accelerate

# Parallel build support (optional)
MPIF90 = mpif90
MPILIBS = 

# Installation (optional, used by 'make install')
PREFIX = /usr/local
BINDIR = $(PREFIX)/bin

# Utility programs
TAR = tar

