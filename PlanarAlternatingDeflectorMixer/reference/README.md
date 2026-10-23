# Reference data — the archived research campaign

Read-only snapshot of the corrected research-fidelity campaign
(`corrected_boundary_v3`, completed 2026-07-26). **No workflow writes here.**
Everything the workflows produce goes to `../results/`, which is never tracked,
so `snakemake ... clean` (`rm -rf results/`) cannot touch these files.

| file | what |
|---|---|
| `corrected_boundary_v3/all_samples.csv` | the twelve screening designs: BO coordinates, CAD parameters, objectives. Read by `docs/exercises/surrogate.py` (assignment tasks 5–6). |
| `corrected_boundary_v3/screening_gate.json` | the screening decision: 12/12 successful, best mixing index 0.1698 against the 0.60 gate — **NO-GO**. |
| `corrected_boundary_v3_baselines/baseline_summary.csv` | the three matched baselines (straight, symmetric deflector, strong alternating) and the analytic straight-channel check. |

These were computed at **research fidelity** (157 604-cell mesh, tight
tolerances). Results produced at workshop fidelity are not comparable and must
never be pooled with them.

A live research campaign (`--campaign research`) writes its own copies of these
files under `../results/`; it does not read or update this folder.
