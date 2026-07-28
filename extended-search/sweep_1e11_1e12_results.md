# Sweep results: 10^11 -> 1.212x10^12

**Run 2026-07-25 to 2026-07-28. Engine `session22_pollard_brent.py`, validated in PR #33/#37.
Final checkpoint `results/session23_checkpoint.json`. Stopped deliberately at 1.212x10^12,
short of the original 10^13 target - see "Why it stopped" below.**

## Two new records

| p | A(p) | d | vs previous |
|---|---|---|---|
| 235,229,251,009 | **367** | 201,394,907 | beats 359 (held since 10^10) |
| 418,383,886,321 | **479** | 2,404,875 | beats 367 |

Both independently re-derived from scratch (not read off the checkpoint): primality
verified, least gateway re-computed by an independent divisor-residue search, divisor
confirmed against `N_A^2`. Full record sequence: 107 -> 167 -> 239 -> 359 -> 367 -> 479.

**Relative to Phase 3a's two candidate models**, this is a genuine surprise in rate, not
in kind. Both models predicted no record at all by 10^13; getting two before even
reaching 1.3x10^12 lands well earlier than either curve implied. Neither is refuted - 479
is comfortably under the saturating model's ceiling of 643 - but the *cadence* of new
records is itself informative and should feed into any future refit of Phase 3a, which
this run does not attempt.

## What this changes about Phase 3a

Nothing conclusively, but it's a real update to the record sequence in
`phase3a_growth_law.md`:

| p | A | gap vs previous |
|---|---|---|
| 3,830,401 | 107 | - |
| 5,462,209 | 167 | 1.4x |
| 32,349,601 | 239 | 5.9x |
| 3,807,728,761 | 359 | 117.7x |
| 235,229,251,009 | 367 | 61.8x |
| 418,383,886,321 | 479 | 1.8x |

The 359->367 gap (61.8x) is in line with the prior escalating pattern; the 367->479 gap
(1.8x) breaks it sharply. Four points was already too few to extrapolate a gap sequence;
six is still too few, but the new data doesn't support "gaps keep growing", which was the
informal expectation stated in the original sweep scoping (`phase3a_sweep_scope.md`).

## Why it stopped short of 10^13

The original scope estimate (PR #32) for the 10^12 -> 10^13 leg was wrong - a units bug
that divided by the worker count twice, understating the true cost by roughly 12x. The
corrected, measured rate on this leg (~14.6 s/segment at `ES_WORKERS=10`) implies the
full run to 10^13 would take **~15 days**, not the ~18-29 hours originally quoted.

Given the original scope's own finding - extending to 10^13 only widens Phase 3a's
`log log X` regressor range by 5%, and is explicitly "low value" for discriminating the
two growth models - 15 days of machine time was not a good trade against two already-found
records. Stopped deliberately at 1.212x10^12 rather than continue un-budgeted.

## State

- `results/session23_checkpoint.json`: `next_lo = 1,212,100,000,001`, `mxA = 479`,
  19,540 witnesses recorded, milestones logged through 10^12.
- Fully resumable: `ES_WORKERS=10 python session22_pollard_brent.py <LIMIT>
  results/session23_checkpoint.json` continues from here if ever revisited.
- No code or engine changes needed to resume - the Miller-Rabin witness set
  `{2,3,5,7,11,13}` remains deterministic up to `N_A ~ 2.5x10^12`, which covers any
  `LIMIT` up to `~1.39x10^13`, per `session22_pollard_brent.py`'s own bound.

## Incidents during this run (both already fixed, PR #33/#37)

1. Sweep initially launched on `session20_corrected_10B.py` (full trial division), cost
   estimated from `session22_pollard_brent.py`'s benchmark - wrong script running,
   ~2.4x slower than expected. Diagnosed live (GPU idle, CPU saturated - not a resource
   problem), switched to the validated fast engine without losing checkpoint state.
2. A local machine crash mid-run; checkpoint and progress survived intact, resumed
   cleanly with no lost work.
3. This document's own cost correction (units bug in the original 10^12->10^13 estimate).

None of the three affected correctness - only cost estimates and wall-clock time.
