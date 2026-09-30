# Git revision stamp for output-file provenance. Generates the source that
# defines frehg::gitRevision() (src/core/GitRevision.hpp), which the run
# record's provenance.git_sha, the HDF5 /frehg2 git_sha attribute and the log
# header read.
#
# Sets FREHG_GIT_REVISION_SOURCE (the generated source: list it in the
# consuming target's sources) and defines the frehg_git_revision target (add
# it as a dependency of that target).
#
# The revision is captured on every build, not at configure time (v2 plan
# V2-A20 finding 5): a commit or a working-tree edit touches no CMake input,
# so a configure-time stamp survived every incremental rebuild and named a
# tree the binary was not built from. The frehg_git_revision target always
# runs FrehgGitRevisionStamp.cmake, which rewrites the source only when the
# stamp changes, so a no-op build recompiles nothing. The stamp is "<sha12>",
# "<sha12>-dirty" or "unknown" (formats: FrehgGitRevisionStamp.cmake); it is
# always the commit SHA, never a tag name. The capture happens as the build
# starts, so an edit made while a build runs shows up in the next build's
# stamp. unit.git_revision_stamp gates all of this.

find_package(Git QUIET)

set(FREHG_GIT_REVISION_SOURCE "${PROJECT_BINARY_DIR}/generated/GitRevision.cpp")
set(_frehg_git_revision_command
  "${CMAKE_COMMAND}"
  "-DFREHG_GIT_EXECUTABLE=${GIT_EXECUTABLE}"
  "-DFREHG_GIT_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
  "-DFREHG_GIT_REVISION_OUTPUT=${FREHG_GIT_REVISION_SOURCE}"
  -P "${CMAKE_CURRENT_LIST_DIR}/FrehgGitRevisionStamp.cmake")

# Once at configure time, so the source exists before the first build (a
# target built without its dependencies, e.g. `make frehg_core/fast`, still
# compiles)...
execute_process(COMMAND ${_frehg_git_revision_command}
  RESULT_VARIABLE _frehg_git_revision_rc)
if(NOT _frehg_git_revision_rc EQUAL 0)
  message(FATAL_ERROR "Frehg2: writing ${FREHG_GIT_REVISION_SOURCE} failed")
endif()
# ...and on every build. BYPRODUCTS lets Ninja re-check the file's timestamp
# after the command runs (restat), so an unchanged stamp recompiles nothing
# there either.
add_custom_target(frehg_git_revision
  COMMAND ${_frehg_git_revision_command}
  BYPRODUCTS "${FREHG_GIT_REVISION_SOURCE}"
  VERBATIM)
