# Sanitizer configuration (FREHG_SANITIZE=ON enables Address + UndefinedBehavior).
#
# Applied globally so that the libraries, the executable, and every test binary
# run under the same instrumentation, as required by the plan's sanitizer lane
# (Section 8.4).

if(FREHG_SANITIZE)
  set(FREHG_SANITIZE_FLAGS -fsanitize=address,undefined -fno-omit-frame-pointer)
  add_compile_options(${FREHG_SANITIZE_FLAGS})
  add_link_options(${FREHG_SANITIZE_FLAGS})
  message(STATUS "Frehg2: Address + UndefinedBehavior sanitizers enabled")
endif()
