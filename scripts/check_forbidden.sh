#!/usr/bin/env bash
# check_forbidden.sh — machine-enforced anti-shortcut gate (upgrade plan §11.3).
#
# Scans the Frehg2 source tree for forbidden patterns and exits nonzero on any
# hit. Run from anywhere; the script locates the repository root itself.

set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="${ROOT}/src"
STATUS=0

fail() {
  echo "FORBIDDEN [$1]:"
  echo "$2"
  STATUS=1
}

scan_src() {
  # $1 = rule name, $2 = extended regex, $3 = extra grep flags (optional)
  local hits
  hits=$(grep -rnE ${3:-} --include='*.cpp' --include='*.hpp' --include='*.h' "$2" "$SRC" 2>/dev/null)
  if [ -n "$hits" ]; then
    fail "$1" "$hits"
  fi
}

scan_src_excluding() {
  # $1 = rule name, $2 = extended regex, $3 = egrep pattern of paths to exempt
  local hits
  hits=$(grep -rnE --include='*.cpp' --include='*.hpp' --include='*.h' "$2" "$SRC" 2>/dev/null | grep -vE "$3")
  if [ -n "$hits" ]; then
    fail "$1" "$hits"
  fi
}

# 1. Deferred-work markers (case-insensitive, comments included).
scan_src "deferred-work markers" \
  'TODO|FIXME|XXX|HACK|WIP|placeholder|stub|not.?implemented|for now|temporar(y|ily)' "-i"

# 2. Disabled code blocks and warning suppression.
scan_src "#if 0 blocks" '^[[:space:]]*#[[:space:]]*if[[:space:]]+0'
scan_src "warning suppression flags" '\-Wno\-'
scan_src "diagnostic pragmas" '#[[:space:]]*pragma[[:space:]]+(GCC|clang)[[:space:]]+diagnostic[[:space:]]+ignored'

# 3. Build-system rules: no globbing, no git submodules.
GLOB_HITS=$(grep -rnE 'file[[:space:]]*\([[:space:]]*GLOB' \
  "$ROOT/CMakeLists.txt" "$ROOT/src" "$ROOT/tests" "$ROOT/cmake" 2>/dev/null | grep -E 'CMakeLists\.txt|\.cmake')
if [ -n "$GLOB_HITS" ]; then
  fail "file(GLOB in CMake" "$GLOB_HITS"
fi
if [ -e "$ROOT/.gitmodules" ]; then
  fail "git submodules" "$ROOT/.gitmodules exists"
fi

# 4. Header hygiene and manual memory management. Deleted special member
#    functions ("= delete") are C++ API design, not memory management, and
#    are exempted.
HDR_NS=$(grep -rnE --include='*.hpp' --include='*.h' 'using[[:space:]]+namespace[[:space:]]+std' "$SRC" 2>/dev/null)
if [ -n "$HDR_NS" ]; then
  fail "using namespace std in headers" "$HDR_NS"
fi
MEM_HITS=$(grep -rnE --include='*.cpp' --include='*.hpp' --include='*.h' \
  '\bnew\b|\bdelete\b|\bmalloc\(|\bfree\(' "$SRC" 2>/dev/null | grep -vE '=[[:space:]]*delete')
if [ -n "$MEM_HITS" ]; then
  fail "raw memory management" "$MEM_HITS"
fi

# 5. Output and process-control discipline.
scan_src_excluding "printf/std::cout outside Logger/main" 'printf\(|std::cout|std::cerr' \
  'core/Logger\.cpp|main\.cpp'
scan_src_excluding "exit/abort outside fatal paths" '\bexit\(|\babort\(' \
  'core/Logger\.cpp|core/PetscSession\.cpp'

# 6. Memory-space and backend-isolation rules.
scan_src "managed memory spaces" 'SharedSpace|CudaUVMSpace'
for moddir in swe gw transport coupling; do
  if [ -d "$SRC/$moddir" ]; then
    MOD_HITS=$(grep -rnE "#[[:space:]]*ifdef[[:space:]]+KOKKOS_ENABLE_(CUDA|HIP)" "$SRC/$moddir" 2>/dev/null)
    if [ -n "$MOD_HITS" ]; then
      fail "backend ifdef in physics module" "$MOD_HITS"
    fi
  fi
done

# 7. Exception-swallowing.
scan_src "catch-all handlers" 'catch[[:space:]]*\([[:space:]]*\.\.\.[[:space:]]*\)'
EMPTY_CATCH=$(grep -rnEz --include='*.cpp' --include='*.hpp' 'catch[[:space:]]*\([^)]*\)[[:space:]]*\{[[:space:]]*\}' "$SRC" 2>/dev/null | tr '\0' '\n')
if [ -n "$EMPTY_CATCH" ]; then
  fail "empty catch bodies" "$EMPTY_CATCH"
fi

# 8. Dropped-feature tripwire (plan §3.2): legacy dead features must not leak in.
scan_src "dropped-feature keywords" \
  'difuwave|diffusive.?wave|subgrid|use_mvg|modified.?picard|newton_iter|waterfall_velocity|pseudo_seepage|evap_model|bctype_SW' "-i"

if [ "$STATUS" -ne 0 ]; then
  echo "check_forbidden.sh: FAILED"
else
  echo "check_forbidden.sh: clean"
fi
exit "$STATUS"
