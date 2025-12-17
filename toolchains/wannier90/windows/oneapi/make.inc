# Wannier90 v3.1.0 make.inc for Windows (Intel oneAPI ifx + MKL)
# Based on Linux ifx template, adapted for Windows
# Uses Intel Fortran (ifx) and Intel MKL

#=====================================================
# For Windows x64 with Intel ifx
#=====================================================
F90=ifx
COMMS=
MPIF90=

# Compiler flags
# -O3: optimization
# -fpp: enable preprocessing
# -Qmkl:sequential: link with Intel MKL sequential libraries
# -warn: enable all warnings (compiler only, not for linker)
FCOPTS=/O3 /fpp /Qmkl:sequential -warn 
# Linker flags (exclude compiler-only flags like -warn)
LDOPTS=/O3 /Qmkl:sequential

#========================================================
# Intel MKL libraries for Windows
# MKLROOT should be set by oneAPI setvars.bat or passed from build script
# For Windows, MKL libraries are in MKLROOT/lib/intel64/
#========================================================


# MKL library directory (use forward slashes for Make)
# On Windows, MKL libraries are directly in lib/, not lib/intel64
MKL_LIB_DIR = $(MKLROOT)/lib

# MKL libraries for sequential build (no OpenMP)
# Interface libraries: mkl_intel_lp64 (64-bit integers, Intel interface)
# Core library: mkl_core
# Threading: mkl_sequential (sequential, no threading)
# Note: On Windows, ifx/ifort will find .lib files automatically with -l flag
LIBS = 
# f90wrap location (optional, for Python interface)
F90WRAP = f90wrap

# Installation (optional, used by 'make install')
PREFIX = C:\wannier90
BINDIR = $(PREFIX)\bin

# Utility programs
TAR = tar

