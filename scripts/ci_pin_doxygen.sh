#!/usr/bin/env bash
# ci_pin_doxygen.sh — install the doxygen version the documentation gate is
# validated against (1.17.0), replacing the runner image's apt doxygen.
#
# ci_build_and_test.sh runs `doxygen docs/Doxyfile` with
# WARN_AS_ERROR = FAIL_ON_WARNINGS. Which constructs count as "documented"
# shifts between doxygen releases: Ubuntu 24.04's apt doxygen is 1.9.8
# (frozen early 2023), and it additionally reports members that are
# documented inside `\name ///@{ ... ///@}` member groups (e.g.
# Grid::comm(), ScalarSolver::surfaceScalar()) as undocumented — false
# positives under the 1.17 parser the project develops against. The gate's
# contract is "the documentation the project validated builds warning-free",
# so CI must run the validated toolchain, not whatever the image froze
# (same reasoning as scripts/ci_fix_mpich_noble.sh for the broken apt
# MPICH). Official release binary from the doxygen GitHub releases,
# sha256-pinned.
#
# Self-gating: a no-op on non-Linux hosts (dev machines get doxygen from
# their package manager; the dev record is brew doxygen 1.17.0) and when the
# PATH doxygen is already the pinned version.

set -euo pipefail

VERSION="1.17.0"
SHA256="75419ef4f446fc1c24ef12514b574e66e898ee6f527c6ae2ad84f91a905823c2"
URL="https://github.com/doxygen/doxygen/releases/download/Release_${VERSION//./_}/doxygen-${VERSION}.linux.bin.tar.gz"

if [ "$(uname -s)" != "Linux" ]; then
  echo "ci_pin_doxygen.sh: not Linux; nothing to do"
  exit 0
fi
current="$(doxygen --version 2>/dev/null | cut -d' ' -f1 || true)"
if [ "$current" = "$VERSION" ]; then
  echo "ci_pin_doxygen.sh: doxygen $current already pinned"
  exit 0
fi
echo "ci_pin_doxygen.sh: replacing doxygen '${current:-none}' with $VERSION"

tmp="$(mktemp -d)"
curl -fsSL --retry 8 --retry-all-errors --retry-delay 5 \
  -o "$tmp/doxygen.tar.gz" "$URL"
echo "$SHA256  $tmp/doxygen.tar.gz" | sha256sum -c -
tar -xzf "$tmp/doxygen.tar.gz" -C "$tmp" "doxygen-${VERSION}/bin/doxygen"
# /usr/local/bin precedes /usr/bin on the runner PATH, so this shadows the
# apt doxygen without uninstalling it.
sudo install -m 0755 "$tmp/doxygen-${VERSION}/bin/doxygen" /usr/local/bin/doxygen
rm -rf "$tmp"

installed="$(doxygen --version | cut -d' ' -f1)"
if [ "$installed" != "$VERSION" ]; then
  echo "ci_pin_doxygen.sh: FATAL — PATH doxygen reports '$installed' after" >&2
  echo "  installing $VERSION to /usr/local/bin (PATH order changed?)." >&2
  exit 1
fi
echo "ci_pin_doxygen.sh: doxygen $installed pinned"
