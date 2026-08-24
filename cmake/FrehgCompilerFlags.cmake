# Strict warning configuration shared by every Frehg2 target.
#
# The flag set below is the zero-warning gate required by the upgrade plan
# (Section 10, P0 exit criteria). Third-party headers are consumed through
# imported targets, whose include directories CMake treats as SYSTEM, so the
# strictness applies to Frehg2 sources only.

add_library(frehg_compiler_flags INTERFACE)
add_library(frehg::compiler_flags ALIAS frehg_compiler_flags)

target_compile_options(frehg_compiler_flags INTERFACE
  -Wall
  -Wextra
  -Wpedantic
  -Wshadow
  -Wconversion)

if(FREHG_WERROR)
  target_compile_options(frehg_compiler_flags INTERFACE -Werror)
endif()
