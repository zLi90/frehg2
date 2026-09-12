#!/usr/bin/env bash
# gpu_acceptance.sh — the p6 GPU acceptance bundle (v2 plan §2B.3).
#
# Run this ON A GPU MACHINE after building frehg2 twice from the same
# commit: once with a CUDA/HIP Kokkos + Kokkos-enabled PETSc (the device
# binary) and once CPU-only (the host reference). It executes the checks CI
# cannot (no GPU on the CI runners) and prints PASS/FAIL per item against
# written criteria, so the result needs no interpretation. A passing bundle
# is what promotes the GPU lanes from 'experimental' to supported
# (v2 plan §2B.4): send back gpu-acceptance-report.txt and the run-record
# files it collects.
#
# Usage:
#   scripts/gpu_acceptance.sh --frehg-gpu PATH --frehg-cpu PATH \
#       --mpiexec PATH [--gpus N] [--work DIR]
#
# Items (criteria in docs/developer-guide/gpu-acceptance.md):
#   1 backend    the device binary really runs on the device backend
#   2 invariance single-GPU result matches 1-rank CPU (volume to 1e-8,
#                gw iterations within 10 %)
#   3 timing     single-GPU vs 1-rank-CPU simulation time (recorded)
#   4 multi-gpu  2 (and --gpus if > 2) ranks, gpu_aware_mpi on and off:
#                volumes match 1-GPU to 1e-8 both ways
#
# Dry-run on a CPU-only machine (CI rehearsal): pass the same binary for
# --frehg-gpu and --frehg-cpu and add --allow-host; item 1 is then reported
# SKIP (host backend) and everything else executes the identical code path.

set -uo pipefail

FREHG_GPU="" FREHG_CPU="" MPIEXEC="" GPUS=2 WORK="gpu-acceptance" ALLOW_HOST=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --frehg-gpu) FREHG_GPU="$2"; shift 2 ;;
    --frehg-cpu) FREHG_CPU="$2"; shift 2 ;;
    --mpiexec)   MPIEXEC="$2"; shift 2 ;;
    --gpus)      GPUS="$2"; shift 2 ;;
    --work)      WORK="$2"; shift 2 ;;
    --allow-host) ALLOW_HOST=1; shift ;;
    *) echo "unknown argument: $1"; exit 2 ;;
  esac
done
[[ -x "$FREHG_GPU" && -x "$FREHG_CPU" && -n "$MPIEXEC" ]] || {
  echo "usage: $0 --frehg-gpu PATH --frehg-cpu PATH --mpiexec PATH [--gpus N] [--work DIR]"
  exit 2
}
# run() changes into each case directory, so binaries must be absolute.
FREHG_GPU="$(cd "$(dirname "$FREHG_GPU")" && pwd)/$(basename "$FREHG_GPU")"
FREHG_CPU="$(cd "$(dirname "$FREHG_CPU")" && pwd)/$(basename "$FREHG_CPU")"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$WORK"
WORK="$(cd "$WORK" && pwd)"
REPORT="$WORK/gpu-acceptance-report.txt"
: > "$REPORT"
PASS_ALL=1
export FI_PROVIDER="${FI_PROVIDER:-tcp}"

# The backend-aware AMG defaults differ host vs device by design (HMIS vs
# PMIS coarsening, v2 plan §2B.2 B2), which would legitimately shift
# iteration counts between the two binaries. For the invariance item both
# runs pin the device-capable configuration on the PETSc command line
# (setOptionDefault yields to it), so any residual drift is a real defect.
PIN_AMG=()
for p in fs gw; do
  PIN_AMG+=("-${p}_pc_hypre_boomeramg_coarsen_type" PMIS
            "-${p}_pc_hypre_boomeramg_relax_type_all" l1scaled-Jacobi)
done

say()  { echo "$*" | tee -a "$REPORT"; }
item() { # item <name> <PASS|FAIL|SKIP|INFO> <detail>
  say "[p6.$1] $2: $3"
  [[ "$2" == "FAIL" ]] && PASS_ALL=0
}

# Stage the >= 1e6-cell synthetic case (the s3/p2 grid) via the harness.
stage() { # stage <tag> -> prints config path
  python3 - "$ROOT" "$WORK" "$1" <<'EOF'
import sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.argv[1]) / "scripts"))
from run_scaling import stage_synthetic
print(stage_synthetic(Path(sys.argv[2]), t_end=60.0, dt=2.0,
                      output_interval=60.0, solver="amg", tag=sys.argv[3]))
EOF
}

run() { # run <binary> <ranks> <config> <log> [extra frehg args...]
  local bin="$1" ranks="$2" config="$3" log="$4"; shift 4
  ( cd "$(dirname "$config")" &&
    OMP_NUM_THREADS=1 "$MPIEXEC" -n "$ranks" "$bin" "$(basename "$config")" "$@" \
      > "$log" 2>&1 )
}

metric() { # metric <case-dir> <field: volume|iters|backend|sim_s>
  python3 - "$1" "$2" <<'EOF'
import sys, yaml
from pathlib import Path
rec = yaml.safe_load((Path(sys.argv[1]) / "out" / "run-record.yaml").read_text())
key = sys.argv[2]
if key == "volume":
    import h5py
    # Surface + subsurface volume: the subsurface term (~1e6 m^3 on this
    # case) keeps the invariance check non-vacuous even while the surface
    # holds only thin films.
    with h5py.File(Path(sys.argv[1]) / "out" / "output.h5") as h:
        surf = float(h["/monitor/mass_audit"][-1, 1])
        sub = float(h["/monitor/gw_mass_audit"][-1, 1])
    print(repr(surf + sub))
elif key == "iters":
    print(rec["solver"]["gw"]["iters_mean"])
elif key == "backend":
    print(rec["provenance"]["kokkos_backend"])
elif key == "sim_s":
    print(rec["timers"]["simulation"]["max_s"])
EOF
}

rel_ok() { # rel_ok <a> <b> <bound> -> exit 0 if |a-b|/max(|a|,eps) <= bound
  python3 - "$1" "$2" "$3" <<'EOF'
import sys
a, b, bound = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
sys.exit(0 if abs(a - b) / max(abs(a), 1e-30) <= bound else 1)
EOF
}

say "p6 GPU acceptance bundle — $(date -u +%Y-%m-%dT%H:%M:%SZ)"
say "gpu binary: $FREHG_GPU"
say "cpu binary: $FREHG_CPU"

# --- item 1: backend -------------------------------------------------------
CFG_GPU="$(stage p6-gpu | tail -1)"
run "$FREHG_GPU" 1 "$CFG_GPU" "$WORK/run-gpu-1.log" "${PIN_AMG[@]}" \
  || { item backend FAIL "single-GPU run exited nonzero (see $WORK/run-gpu-1.log)"; }
BACKEND="$(metric "$(dirname "$CFG_GPU")" backend || echo unknown)"
case "$BACKEND" in
  *Cuda*|*HIP*|*SYCL*) item backend PASS "device backend '$BACKEND'" ;;
  *) if [[ "$ALLOW_HOST" == 1 ]]; then
       item backend SKIP "host backend '$BACKEND' (--allow-host dry run)"
     else
       item backend FAIL "expected a device backend, run record says '$BACKEND'"
     fi ;;
esac
V_GPU="$(metric "$(dirname "$CFG_GPU")" volume)"
I_GPU="$(metric "$(dirname "$CFG_GPU")" iters)"
T_GPU="$(metric "$(dirname "$CFG_GPU")" sim_s)"
cp "$(dirname "$CFG_GPU")/out/run-record.yaml" "$WORK/record-gpu-1.yaml"

# --- item 2 + 3: invariance and timing vs 1-rank CPU -----------------------
CFG_CPU="$(stage p6-cpu | tail -1)"
run "$FREHG_CPU" 1 "$CFG_CPU" "$WORK/run-cpu-1.log" "${PIN_AMG[@]}" \
  || item invariance FAIL "CPU reference run exited nonzero"
V_CPU="$(metric "$(dirname "$CFG_CPU")" volume)"
I_CPU="$(metric "$(dirname "$CFG_CPU")" iters)"
T_CPU="$(metric "$(dirname "$CFG_CPU")" sim_s)"
cp "$(dirname "$CFG_CPU")/out/run-record.yaml" "$WORK/record-cpu-1.yaml"
if rel_ok "$V_CPU" "$V_GPU" 1e-8; then
  item invariance PASS "final volume GPU vs CPU rel diff <= 1e-8 ($V_GPU vs $V_CPU)"
else
  item invariance FAIL "final volume GPU $V_GPU vs CPU $V_CPU exceeds 1e-8 relative"
fi
if rel_ok "$I_CPU" "$I_GPU" 0.10; then
  item invariance PASS "gw iterations within 10 % (GPU $I_GPU, CPU $I_CPU)"
else
  item invariance FAIL "gw iterations drift > 10 % (GPU $I_GPU, CPU $I_CPU)"
fi
item timing INFO "simulation time: GPU $T_GPU s vs CPU(1 rank) $T_CPU s"

# --- item 4: multi-GPU, gpu_aware_mpi on and off ---------------------------
for mode in on off; do
  for n in 2 $([[ "$GPUS" -gt 2 ]] && echo "$GPUS"); do
    CFG_N="$(stage "p6-n${n}-${mode}" | tail -1)"
    python3 - "$CFG_N" "$mode" <<'EOF'
import sys, yaml
doc = yaml.safe_load(open(sys.argv[1]))
doc["runtime"] = {"gpu_aware_mpi": sys.argv[2]}
yaml.safe_dump(doc, open(sys.argv[1], "w"), sort_keys=False)
EOF
    if run "$FREHG_GPU" "$n" "$CFG_N" "$WORK/run-gpu-${n}-${mode}.log" "${PIN_AMG[@]}"; then
      V_N="$(metric "$(dirname "$CFG_N")" volume)"
      T_N="$(metric "$(dirname "$CFG_N")" sim_s)"
      if rel_ok "$V_GPU" "$V_N" 1e-8; then
        item multi-gpu PASS "n=$n gpu_aware_mpi=$mode: volume matches 1-GPU (sim ${T_N}s)"
      else
        item multi-gpu FAIL "n=$n gpu_aware_mpi=$mode: volume $V_N vs 1-GPU $V_GPU exceeds 1e-8"
      fi
      cp "$(dirname "$CFG_N")/out/run-record.yaml" "$WORK/record-gpu-${n}-${mode}.yaml"
    else
      item multi-gpu FAIL "n=$n gpu_aware_mpi=$mode exited nonzero (see log)"
    fi
  done
done

say ""
if [[ "$PASS_ALL" == 1 ]]; then
  say "p6 BUNDLE: PASS — send $REPORT and $WORK/record-*.yaml back; the GPU"
  say "lanes are promoted from 'experimental' on receipt (v2 plan §2B.4)."
else
  say "p6 BUNDLE: FAIL — see the FAIL items above; include $WORK/*.log with the report."
fi
[[ "$PASS_ALL" == 1 ]]
