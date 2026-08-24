# Frehg2 benchmark configurations

One directory per validation benchmark (plan §2.3, §9). Each holds the
frozen-v2 YAML configuration(s), the input data converted once from
`legacy/benchmarks/` (provenance in the YAML headers and per-case READMEs),
and — for b1/b2 — the archival v1 YAML drafts consumed by
`tools/migrate_yaml_v1_to_v2.py`.

| Case | Config(s) | Modules | Gate phase |
|---|---|---|---|
| b1-sw | `b1-sw/b1-sw.yaml` | SWE | P1 |
| b2-gw | `b2-gw/b2-gw.yaml` | GW | P2 |
| b3-kirkland | `b3-kirkland/b3-kirkland.yaml` | GW | P2 |
| b4-govindaraju | `b4-govindaraju/b4-govindaraju.yaml` | SWE (Chezy) | P1 |
| b5-vcatchment | `b5-vcatchment/b5-vcatchment{,-norain}.yaml` | SWE+GW | P3 |
| b6-kuan | `b6-kuan/b6-kuan-{ss,td}.yaml` | SWE+GW+transport | P4 |

Golden reference data is **not** committed here; the regression harness
(P1+) reads it from a `legacy/benchmarks/` checkout — a development-side
archive (legacy model outputs, digitized reference curves, plotting
scripts) that is **not distributed with this repository** (plan §8.3).
Without it the golden-comparison regression tests cannot run; everything
else (the unit and mpi test labels, `--validate`, and running every
benchmark and validation case) is fully self-contained.

All configurations are validated by CI on every commit
(`ctest -R validate.`), and `frehg --validate <case>.yaml` reproduces that
check locally.

Input data conventions (see `src/io/GridDataReader.hpp`):
- rasters: headered (`ncols/nrows/...`) or flat lists; the first data row is
  j = 0 (file order — the convention the legacy references were consumed
  with, not ESRI north-up);
- 3D fields (soil ids, subsurface ICs): flat lists ordered `(j*nx+i)*nz+k`;
- time series: two-column `t value` lines, seconds and SI units.
