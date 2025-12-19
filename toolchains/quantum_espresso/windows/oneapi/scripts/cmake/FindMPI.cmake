# FindMPI.cmake shim for MS-MPI on Windows
# This shim bypasses CMake's built-in FindMPI.cmake to avoid MPI_Fortran_WORKS try_compile tests
# that fail with Intel ifx + MS-MPI on Windows.

# Locate MS-MPI SDK
# Allow override via MSMPI_SDK environment variable, otherwise use default location
if (DEFINED ENV{MSMPI_SDK})
    set(MSMPI_SDK_ROOT "$ENV{MSMPI_SDK}")
else()
    set(MSMPI_SDK_ROOT "C:/Program Files (x86)/Microsoft SDKs/MPI")
endif()

# Convert to CMake path format (forward slashes)
file(TO_CMAKE_PATH "${MSMPI_SDK_ROOT}" MSMPI_SDK_ROOT)

# Set paths
set(MSMPI_INCLUDE_DIR_BASE "${MSMPI_SDK_ROOT}/Include")
set(MSMPI_INCLUDE_DIR_ARCH "${MSMPI_SDK_ROOT}/Include/x64")
set(MSMPI_INCLUDE_DIRS "${MSMPI_INCLUDE_DIR_BASE};${MSMPI_INCLUDE_DIR_ARCH}")
set(MSMPI_LIB_DIR "${MSMPI_SDK_ROOT}/Lib/x64")
set(MSMPI_LIB "${MSMPI_LIB_DIR}/msmpi.lib")
set(MSMPI_FEC_LIB "${MSMPI_LIB_DIR}/msmpifec.lib")

# Locate mpiexec
set(MSMPI_MPIEXEC "C:/Program Files/Microsoft MPI/Bin/mpiexec.exe")
file(TO_CMAKE_PATH "${MSMPI_MPIEXEC}" MSMPI_MPIEXEC)

# Validate required files exist
if (NOT EXISTS "${MSMPI_INCLUDE_DIR_BASE}/mpif.h")
    message(FATAL_ERROR "MS-MPI shim: mpif.h not found at ${MSMPI_INCLUDE_DIR_BASE}/mpif.h")
endif()

if (NOT EXISTS "${MSMPI_LIB}")
    message(FATAL_ERROR "MS-MPI shim: msmpi.lib not found at ${MSMPI_LIB}. Library directory: ${MSMPI_LIB_DIR}")
endif()

if (NOT EXISTS "${MSMPI_FEC_LIB}")
    message(FATAL_ERROR "MS-MPI shim: msmpifec.lib not found at ${MSMPI_FEC_LIB}. Library directory: ${MSMPI_LIB_DIR}")
endif()

if (NOT EXISTS "${MSMPI_MPIEXEC}")
    message(FATAL_ERROR "MS-MPI shim: mpiexec.exe not found at ${MSMPI_MPIEXEC}")
endif()

# Set standard FindMPI outputs
set(MPI_FOUND TRUE CACHE BOOL "MPI found" FORCE)
set(MPI_C_FOUND TRUE CACHE BOOL "MPI C found" FORCE)
set(MPI_Fortran_FOUND TRUE CACHE BOOL "MPI Fortran found" FORCE)
set(MPI_Fortran_HAVE_F90_MODULE FALSE)
set(MPI_Fortran_HAVE_F77_HEADER TRUE)
set(MPI_C_INCLUDE_DIRS "${MSMPI_INCLUDE_DIRS}" CACHE STRING "MPI C include directories" FORCE)
set(MPI_Fortran_INCLUDE_DIRS "${MSMPI_INCLUDE_DIRS}" CACHE STRING "MPI Fortran include directories" FORCE)
set(MPI_C_LIBRARIES "${MSMPI_LIB}" CACHE STRING "MPI C libraries" FORCE)
set(MPI_Fortran_LIBRARIES "${MSMPI_FEC_LIB};${MSMPI_LIB}" CACHE STRING "MPI Fortran libraries" FORCE)
set(MPIEXEC_EXECUTABLE "${MSMPI_MPIEXEC}" CACHE FILEPATH "MPI executable" FORCE)
set(MPIEXEC_NUMPROC_FLAG "-n")
set(MPIEXEC_PREFLAGS "")

# Create/update imported targets expected by QE
# MPI::MPI_Fortran - must link msmpifec.lib and msmpi.lib
if (NOT TARGET MPI::MPI_Fortran)
    add_library(MPI::MPI_Fortran INTERFACE IMPORTED)
endif()
set_target_properties(MPI::MPI_Fortran PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${MSMPI_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${MSMPI_FEC_LIB};${MSMPI_LIB}"
)

# MPI::MPI_C - optional but good hygiene
if (NOT TARGET MPI::MPI_C)
    add_library(MPI::MPI_C INTERFACE IMPORTED)
endif()
set_target_properties(MPI::MPI_C PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${MSMPI_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${MSMPI_LIB}"
)

# Print status
message(STATUS "Using MS-MPI shim: includes=${MSMPI_INCLUDE_DIRS} libs=msmpifec;msmpi")

