#!/usr/bin/env bash
# ci_fix_mpich_noble.sh — replace Ubuntu 24.04 (noble)'s broken MPICH with the
# fixed 4.2.0-5.1 build, then prove multi-rank launch actually works.
#
# Noble ships mpich 4.2.0-5build3 built against PMIx, but its bundled Hydra
# mpiexec speaks only classic PMI-1 — launched processes cannot obtain their
# process-group identity, so EVERY rank silently degrades to a size-1
# MPI_COMM_WORLD ("mpiexec -n 4" runs four independent singletons). For frehg
# that meant the mpi ctest label never really ran multi-rank in CI: the halo
# drivers passed trivially on their single-rank paths, and the parallel-HDF5
# driver's N singletons raced on one shared filename and flaked
# (mpi.core.n4). Upstream references:
#   https://bugs.launchpad.net/ubuntu/+source/mpich/+bug/2072338
#   https://github.com/pmodels/mpich/issues/7064
#   https://github.com/actions/runner-images/issues/13204
# No environment-variable workaround exists; the fix is the 4.2.0-5.1 rebuild
# (PMIx linkage dropped, same upstream 4.2.0 / same libmpich.so.12 ABI), which
# never landed in noble-updates — so it is installed here from Launchpad's
# official Ubuntu build farm, sha256-pinned. Because the soname and ABI are
# unchanged, the PETSc/HDF5 stacks built against 5build3 (including the cached
# dependency prefix) keep working untouched — no cache-key bump.
#
# Self-gating: a no-op on non-dpkg systems (macOS dev machines) and on any
# image whose MPICH is already >= 4.2.0-5.1. Run it in every workflow right
# after the apt install step, before anything launches mpiexec.

set -euo pipefail

if ! command -v dpkg >/dev/null 2>&1; then
  echo "ci_fix_mpich_noble.sh: not a dpkg system; nothing to do"
  exit 0
fi
installed="$(dpkg-query -W -f '${Version}' libmpich12 2>/dev/null || true)"
if [ -z "$installed" ]; then
  echo "ci_fix_mpich_noble.sh: libmpich12 not installed; nothing to do"
  exit 0
fi

FIXED="4.2.0-5.1"
BASE="https://launchpad.net/ubuntu/+source/mpich/4.2.0-5.1/+build/28285882/+files"

if dpkg --compare-versions "$installed" ge "$FIXED"; then
  echo "ci_fix_mpich_noble.sh: libmpich12 $installed >= $FIXED; package already fixed"
else
  echo "ci_fix_mpich_noble.sh: libmpich12 $installed is the broken PMIx build; installing $FIXED"
  tmp="$(mktemp -d)"
  # sha256 of the official Launchpad amd64 builds (build 28285882), pinned so a
  # tampered or truncated download cannot install. Launchpad 502s transiently;
  # retry hard.
  fetch() { # fetch <file> <sha256>
    curl -fsSL --retry 8 --retry-all-errors --retry-delay 5 \
      -o "$tmp/$1" "$BASE/$1"
    echo "$2  $tmp/$1" | sha256sum -c -
  }
  fetch "libmpich12_${FIXED}_amd64.deb"   "13a233b7dd6d8a51eccf58d981cdb1ae110b673e4dfdeb607fe439918340238f"
  fetch "mpich_${FIXED}_amd64.deb"        "c696dff0975f1857faaae0e450ce6262d13693a7894e523741e6b7581beee0ad"
  fetch "libmpich-dev_${FIXED}_amd64.deb" "35cebdadbe718a7ee0fa402d0db157e5a471e3d373822b6f42445019be33c4b2"
  # One dpkg transaction: libmpich-dev/mpich depend on libmpich12 (= 4.2.0-5.1)
  # exactly, so the three must go in together. Runtime deps (hwloc, ucx, slurm,
  # gfortran) are already satisfied by the apt install of the 5build3 packages
  # this replaces — identical dependency sets. Do NOT "apt-get install -f" on
  # failure: apt would resolve it by reverting to the broken 5build3.
  sudo dpkg -i \
    "$tmp/libmpich12_${FIXED}_amd64.deb" \
    "$tmp/mpich_${FIXED}_amd64.deb" \
    "$tmp/libmpich-dev_${FIXED}_amd64.deb"
  # Keep a later apt operation from silently downgrading back to 5build3.
  sudo apt-mark hold libmpich12 mpich libmpich-dev >/dev/null 2>&1 || true
  rm -rf "$tmp"
  dpkg-query -W libmpich12 mpich libmpich-dev
fi

# Prove it, whatever path got here: compile a two-rank hello and require that
# the launched world really has size 2. This is the exact failure mode the
# broken package produces, so a regression (image update, cache of debs gone
# stale, a future re-break) fails THIS step with an unambiguous message
# instead of resurfacing as a flaky parallel-HDF5 race an hour later.
# Same runner pins the test environment uses: UCX otherwise probes InfiniBand
# verbs at MPI_Init and aborts on RDMA-less runners.
export UCX_TLS="${UCX_TLS:-tcp,self,sm}"
export FI_PROVIDER="${FI_PROVIDER:-tcp}"
verify="$(mktemp -d)"
cat > "$verify/hello.c" <<'EOF'
#include <mpi.h>
#include <stdio.h>
int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = -1, size = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  printf("rank %d of %d\n", rank, size);
  MPI_Finalize();
  return size == 2 ? 0 : 1;
}
EOF
MPICC="$(command -v mpicc.mpich || command -v mpicc)"
MPIEXEC="$(command -v mpiexec.mpich || command -v mpiexec)"
"$MPICC" -O0 -o "$verify/hello" "$verify/hello.c"
# timeout: a wedged PMI handshake must fail this step, not hang the job.
if timeout 120 "$MPIEXEC" -n 2 "$verify/hello"; then
  echo "ci_fix_mpich_noble.sh: verified — mpiexec -n 2 forms a size-2 MPI_COMM_WORLD"
  rm -rf "$verify"
else
  echo "ci_fix_mpich_noble.sh: FATAL — mpiexec -n 2 still launches singleton" >&2
  echo "  ranks (size-1 MPI_COMM_WORLD each): the PMI handshake between the" >&2
  echo "  hydra launcher and libmpich is broken on this image. Every mpi-label" >&2
  echo "  test would run as racing singletons. Refusing to continue." >&2
  exit 1
fi
