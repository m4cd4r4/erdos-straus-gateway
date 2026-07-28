# Tier 1 verification: every A>=100 witness, independently re-derived

**Run 2026-07-28. Script `calibrate_tier1_verify.py`, exhaustive full-population
run inline (not the calibration sample). Raw calibration output in
`results/calibrate_tier1_results.json`.**

## What this closes

The paper's verification-standard Remark (Discussion, §6.1) previously said the
extended search past `10^11` had "only the two record witnesses --- not the
classification of every intervening prime --- independently re-derived." That
was accurate but conservative: it undersold what's now been done.

**All 19,540 witnesses with `A >= 100` across the full run (`10^6` through
`1.212 x 10^12`) have been independently re-verified, exhaustively, not
sampled.** For each: primality of `p`, that `d` genuinely divides `N_A^2`,
that `d = t_A mod A`, and - the expensive part - that **every smaller
candidate gateway `B < A` fails**, confirmed by full divisor-residue
enumeration (not a partial check) for each one.

Zero failures. Max `A` independently confirmed at `479`, matching the primary
engine exactly.

## Why this was worth doing immediately rather than staying a calibration

The original ask was a cost calibration to decide whether a full Tier-1 pass
was worth running. The calibration sample (100 witnesses, stratified) came in
at ~62 projected core-seconds for the full population - far cheaper than a
back-of-envelope estimate assuming the largest-known record (`A=479`, 46
gateway checks) was representative. It isn't: the population's median `A` is
107, not 479, so the typical witness needs ~14 gateway checks, each a fast
factorisation of an `~10^11`-scale integer plus a residue-set build over its
divisors - milliseconds, not seconds.

Two buckets (`A` in `[250,400)` and `[400,500)`) had no coverage in the
100-witness sample - only 19 witnesses exist there in the whole population,
too few for random sampling to land on reliably. Rather than extrapolate
across a gap that happens to contain both actual records, all 19 were run
directly: 0.14 seconds total, all passing, including independent
re-confirmation of the original `A=359` record (`p=3,807,728,761`) alongside
the two new ones.

Given the full high-`A` tail costs under a minute, running the complete
19,540-witness population was cheaper than continuing to reason about whether
to run it.

## What this does NOT close

This is Tier 1 of three tiers considered (see the scoping discussion this
followed). It confirms the **tail** - every prime whose gateway search was
expensive enough to produce a witness. It does not touch the aggregate stage
counts (`primes`, `thm1`, `thm2`, `thm4`, `caseA`, `nqr7`, `gw` totals) or the
full `A`-distribution for primes solved by small `A` (7, 11, 19, 23, ...),
which is the vast majority of the ~700 million primes in `[10^11, 1.212 x
10^12]`. Matching the `<=10^11` table's actual standard - independent
cross-validation of *every* prime's classification - is Tier 2, and nothing
here estimates its cost: Tier 1's cheapness doesn't generalise, since
verifying an easy prime (solved at `A=7`) needs a correct factorisation and a
single divisor check, not the exhaustive multi-gateway search the tail
required, so per-prime cost and total prime count are both completely
different regimes.

## Reproduction

```
cd extended-search
python calibrate_tier1_verify.py 100     # the calibration (fast, ~1s)
python verify_all_witnesses.py           # the full population (~1 min)
```

`results/tier1_full_verification.json` records the outcome (19,540 / 19,540
pass, `max_A_confirmed = 479`) for anyone re-running it.
