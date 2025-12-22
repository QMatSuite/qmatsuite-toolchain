set(GITREV_BARE_FILE git-rev.h)
set(GITREV_FILE ${qe_BINARY_DIR}/${GITREV_BARE_FILE})

# Generate git-rev.h without shell pipes/redirection for cross-platform builds
set(GITREV_SCRIPT ${qe_BINARY_DIR}/cmake/git_rev_generate.cmake)
file(MAKE_DIRECTORY ${qe_BINARY_DIR}/cmake)
set(_git_dollar "$")
file(WRITE ${GITREV_SCRIPT}
    "cmake_minimum_required(VERSION 3.16)\n"
    "set(branch \"UNKNOWN\")\n"
    "set(hash \"UNKNOWN\")\n"
    "set(date \"UNKNOWN\")\n"
    "set(subject \"UNKNOWN\")\n"
    "if(EXISTS \"${_git_dollar}{qe_SOURCE_DIR}/.git\" AND GIT_EXECUTABLE)\n"
    "  execute_process(COMMAND \"${_git_dollar}{GIT_EXECUTABLE}\" rev-parse --abbrev-ref HEAD\n"
    "    OUTPUT_VARIABLE branch\n"
    "    OUTPUT_STRIP_TRAILING_WHITESPACE\n"
    "    RESULT_VARIABLE rv_branch)\n"
    "  if(rv_branch)\n"
    "    set(branch \"UNKNOWN\")\n"
    "  endif()\n"
    "  execute_process(COMMAND \"${_git_dollar}{GIT_EXECUTABLE}\" describe --always --dirty --abbrev=40 --match=NoTagWithThisName\n"
    "    OUTPUT_VARIABLE hash_full\n"
    "    OUTPUT_STRIP_TRAILING_WHITESPACE\n"
    "    RESULT_VARIABLE rv_hash)\n"
    "  if(NOT rv_hash)\n"
    "    set(hash \"${_git_dollar}{hash_full}\")\n"
    "  else()\n"
    "    execute_process(COMMAND \"${_git_dollar}{GIT_EXECUTABLE}\" rev-parse --short HEAD\n"
    "      OUTPUT_VARIABLE hash\n"
    "      OUTPUT_STRIP_TRAILING_WHITESPACE\n"
    "      RESULT_VARIABLE rv_hash_short)\n"
    "    if(rv_hash_short)\n"
    "      set(hash \"UNKNOWN\")\n"
    "    endif()\n"
    "  endif()\n"
    "  execute_process(COMMAND \"${_git_dollar}{GIT_EXECUTABLE}\" log -1 --format=%ad\n"
    "    OUTPUT_VARIABLE date\n"
    "    OUTPUT_STRIP_TRAILING_WHITESPACE\n"
    "    RESULT_VARIABLE rv_date)\n"
    "  if(rv_date)\n"
    "    set(date \"UNKNOWN\")\n"
    "  endif()\n"
    "  execute_process(COMMAND \"${_git_dollar}{GIT_EXECUTABLE}\" log -1 --format=%s\n"
    "    OUTPUT_VARIABLE subject\n"
    "    OUTPUT_STRIP_TRAILING_WHITESPACE\n"
    "    RESULT_VARIABLE rv_subject)\n"
    "  if(rv_subject)\n"
    "    set(subject \"UNKNOWN\")\n"
    "  endif()\n"
    "endif()\n"
    "string(REPLACE \"\\\"\" \"\\\\\\\"\" subject_escaped \"${_git_dollar}{subject}\")\n"
    "file(WRITE \"${_git_dollar}{GITREV_FILE}\"\n"
    "  \"#define GIT_BRANCH_RAW \\\"${_git_dollar}{branch}\\\"\\n\"\n"
    "  \"#define GIT_HASH_RAW \\\"${_git_dollar}{hash}\\\"\\n\"\n"
    "  \"#define GIT_COMMIT_LAST_CHANGED_RAW \\\"${_git_dollar}{date}\\\"\\n\"\n"
    "  \"#define GIT_COMMIT_SUBJECT_RAW \\\"${_git_dollar}{subject_escaped}\\\"\\n\")\n")

add_custom_target(gitrev
    COMMAND ${CMAKE_COMMAND}
        -DGIT_EXECUTABLE=${GIT_EXECUTABLE}
        -Dqe_SOURCE_DIR=${qe_SOURCE_DIR}
        -DGITREV_FILE=${GITREV_FILE}
        -P ${GITREV_SCRIPT}
    WORKING_DIRECTORY ${qe_SOURCE_DIR}
    VERBATIM)

# Print some configure-time git info without shell pipelines
set(_git_branch "UNKNOWN")
set(_git_hash "UNKNOWN")
if(EXISTS "${qe_SOURCE_DIR}/.git" AND GIT_EXECUTABLE)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
        OUTPUT_VARIABLE _git_branch
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _rv_branch
        WORKING_DIRECTORY ${qe_SOURCE_DIR})
    if(_rv_branch)
        set(_git_branch "UNKNOWN")
    endif()
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --always --dirty --abbrev=40 --match="NoTagWithThisName"
        OUTPUT_VARIABLE _git_hash
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _rv_hash
        WORKING_DIRECTORY ${qe_SOURCE_DIR})
    if(_rv_hash)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
            OUTPUT_VARIABLE _git_hash
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE _rv_hash_short
            WORKING_DIRECTORY ${qe_SOURCE_DIR})
        if(_rv_hash_short)
            set(_git_hash "UNKNOWN")
        endif()
    endif()
endif()
message("   Git branch: ${_git_branch}")
message("   Git commit hash: ${_git_hash}")
