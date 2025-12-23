# qe-win-versioninfo.cmake
# Add Windows VERSIONINFO resource to ALL executable targets after configure.
#
# Use:
#   -DCMAKE_PROJECT_TOP_LEVEL_INCLUDES=".../qe-win-versioninfo.cmake"
#   -DQMS_EMBED_VERSIONINFO=ON
#   -DQMS_PRODUCT="Quantum ESPRESSO"
#   -DQMS_COMPANY="QMatSuite"
#   -DQMS_VERSION_STR="7.5.0"
#   -DQMS_BUILD_FLAVOR="Windows oneAPI/MS-MPI"

if(NOT WIN32)
  return()
endif()

if(DEFINED QMS_EMBED_VERSIONINFO AND NOT QMS_EMBED_VERSIONINFO)
  return()
endif()

set(QMS_COMPANY      "QMatSuite" CACHE STRING "CompanyName for VERSIONINFO")
set(QMS_PRODUCT      "Quantum ESPRESSO" CACHE STRING "ProductName for VERSIONINFO")
set(QMS_VERSION_STR  "7.5.0" CACHE STRING "ProductVersion/FileVersion string")
set(QMS_COPYRIGHT    "© 2025" CACHE STRING "Copyright string")
set(QMS_BUILD_FLAVOR "" CACHE STRING "Extra flavor appended to FileDescription")

# Parse QMS_VERSION_STR -> 4 integers (a,b,c,d)
set(_qms_ver "${QMS_VERSION_STR}")
string(REPLACE "." ";" _qms_ver_list "${_qms_ver}")
list(LENGTH _qms_ver_list _qms_ver_len)

set(_qms_v0 0)
set(_qms_v1 0)
set(_qms_v2 0)
set(_qms_v3 0)

if(_qms_ver_len GREATER 0)
  list(GET _qms_ver_list 0 _qms_v0)
endif()
if(_qms_ver_len GREATER 1)
  list(GET _qms_ver_list 1 _qms_v1)
endif()
if(_qms_ver_len GREATER 2)
  list(GET _qms_ver_list 2 _qms_v2)
endif()
if(_qms_ver_len GREATER 3)
  list(GET _qms_ver_list 3 _qms_v3)
endif()

# RC template written once into build dir
set(_qms_rc_in "${CMAKE_BINARY_DIR}/qms_versioninfo.rc.in")
if(NOT EXISTS "${_qms_rc_in}")
  file(WRITE "${_qms_rc_in}" [=[
#include <winver.h>

VS_VERSION_INFO VERSIONINFO
 FILEVERSION     @QMS_V0@,@QMS_V1@,@QMS_V2@,@QMS_V3@
 PRODUCTVERSION  @QMS_V0@,@QMS_V1@,@QMS_V2@,@QMS_V3@
 FILEFLAGSMASK   0x3fL
 FILEFLAGS       0x0L
 FILEOS          0x40004L
 FILETYPE        0x1L
 FILESUBTYPE     0x0L
BEGIN
  BLOCK "StringFileInfo"
  BEGIN
    BLOCK "040904B0"
    BEGIN
      VALUE "CompanyName",      "@QMS_COMPANY@"
      VALUE "ProductName",      "@QMS_PRODUCT@"
      VALUE "FileDescription",  "@QMS_DESCRIPTION@"
      VALUE "FileVersion",      "@QMS_VERSION_STR@"
      VALUE "ProductVersion",   "@QMS_VERSION_STR@"
      VALUE "OriginalFilename", "@QMS_ORIGINAL_FILENAME@"
      VALUE "InternalName",     "@QMS_INTERNAL_NAME@"
      VALUE "LegalCopyright",   "@QMS_COPYRIGHT@"
    END
  END
  BLOCK "VarFileInfo"
  BEGIN
    VALUE "Translation", 0x0409, 1200
  END
END
]=])
endif()

# Function that will run at the end of the top-level configure
function(qms_attach_versioninfo_to_all_exes)
  # Collect all targets in the entire buildsystem
  get_property(_allTargets GLOBAL PROPERTY TARGETS)
  if(NOT _allTargets)
    return()
  endif()

  foreach(_tgt IN LISTS _allTargets)
    if(NOT TARGET "${_tgt}")
      continue()
    endif()

    get_target_property(_type "${_tgt}" TYPE)
    if(NOT _type STREQUAL "EXECUTABLE")
      continue()
    endif()

    # Skip imported executables (rare, but be safe)
    get_target_property(_imported "${_tgt}" IMPORTED)
    if(_imported)
      continue()
    endif()

    # Generate per-target rc
    set(QMS_V0 "${_qms_v0}")
    set(QMS_V1 "${_qms_v1}")
    set(QMS_V2 "${_qms_v2}")
    set(QMS_V3 "${_qms_v3}")

    if(QMS_BUILD_FLAVOR STREQUAL "")
      set(QMS_DESCRIPTION "${_tgt}")
    else()
      set(QMS_DESCRIPTION "${_tgt} (${QMS_BUILD_FLAVOR})")
    endif()

    set(QMS_ORIGINAL_FILENAME "${_tgt}.exe")
    set(QMS_INTERNAL_NAME "${_tgt}")

    set(_rc_out "${CMAKE_BINARY_DIR}/qms_versioninfo_${_tgt}.rc")
    configure_file("${_qms_rc_in}" "${_rc_out}" @ONLY)

    # Attach
    target_sources("${_tgt}" PRIVATE "${_rc_out}")
  endforeach()
endfunction()

# Defer to after project() and the full configure has defined targets (avoids messing with FindGit etc.)
cmake_language(DEFER CALL qms_attach_versioninfo_to_all_exes)
