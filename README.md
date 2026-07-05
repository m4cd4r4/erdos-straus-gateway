# Bounded Gateway Parameters for the Erdos-Straus Conjecture

**Gateway Decompositions via Divisors of N^2**

![Results overview: prime classification, A-value distribution, and max A stabilisation](banner.png)

## Overview

The [Erdos-Straus conjecture](https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93Straus_conjecture) (1948) asserts that for every integer n >= 2, the equation

    4/n = 1/x + 1/y + 1/z

has a solution in positive integers x, y, z.

This repository contains the paper, verification code, and figures for a unified algebraic approach to the prime case. Every prime `p = 3 (mod 4)` is handled classically with `A = 1`. For every prime `p = 1 (mod 4)`, a single auxiliary parameter `A` produces an explicit decomposition: three algebraic propositions cover all but a density-zero residual class, and a short search over prime values of `A` resolves the rest.

Verification to **10^11** (all 4,118,054,813 primes) succeeds with just **32 values of A** for the `p = 1 (mod 4)` primes - the value `A = 3` for the algebraic classes plus 31 distinct gateway values, all `= 3 (mod 4)` and all `<= 359`.

## Key Result

**Theorem.** For every prime `p = 1 (mod 4)` in the residual class, there exists a prime `A <= 359` with `A = 3 (mod 4)`, `4 | (p + A)`, and a divisor `d` of `((p+A)/4)^2` satisfying `d = -p^2 * 4^{-1} (mod A)`, such that

    x = (p+A)/4,  y = (p*x + d)/A,  z = p*x*y/d

are positive integers giving `4/p = 1/x + 1/y + 1/z`.

The critical insight is that the integrality condition requires `d | N^2` (where `N = (p+A)/4`), which is **strictly weaker** than `d | N`. A nontrivial fraction of the hardest primes use a divisor `d > N` that is invisible under the stronger condition.

## Bounded-A Phenomenon

The maximum `A` grows extremely slowly with scale, and does not grow at all across the most recent decade:

| Limit | max A | Increase |
|-------|-------|----------|
| 10^6 | 79 | - |
| 10^7 | 167 | +88 |
| 10^8 | 239 | +72 |
| 10^9 | 239 | +0 |
| 10^10 | 359 | +120 |
| 10^11 | 359 | +0 |

The record-setting prime is `p = 3,807,728,761` (`A = 359`, `d = 1935`), which lies **below 10^10** - nothing in `(10^10, 10^11]` beats it. Across 4.1 billion primes, only **47** require `A >= 199` and only **5** require `A >= 251`. This slow, eventually-flat growth motivates the **bounded-A conjecture**: `max A` is absolutely bounded. If true, the full Erdos-Straus conjecture follows.

## Paper

The main paper (`erdos_straus_gateway.tex`) contains:

- **Section 3**: The Gateway Decomposition theorem and the `d | N^2` lemma
- **Section 4**: Algebraic existence proofs covering ~97.1% of all primes (at 10^6, rising with the bound)
- **Section 5**: Computational verification for the remaining residual class (Case B QR7 primes), including an *Anatomy of the hardest prime* - a discrete-log proof of why `p = 3,807,728,761` forces `A = 359` (its N^2 divisors fill every unit mod 359)
- **Section 6**: Discussion of the bounded-A phenomenon and the path to a full proof

The residual class - `p = 1 (mod 24)`, a quadratic residue mod 7, with every prime factor of `(p+3)/4` congruent to `1 (mod 3)` - is related to but **not identical with** the classical Mordell subset `{1, 121, 169, 289, 361, 529} (mod 840)`; Remark 4.5 distinguishes them with explicit examples.

## Verification

Run the self-contained verification script (Python 3.6+, standard library only):

```bash
python verify.py              # all primes to 10^6  (~10 seconds)
python verify.py 10000000     # to 10^7             (~2 minutes)
python verify.py 100000000    # to 10^8             (~30 minutes)
python verify.py 1000000000   # to 10^9             (~90 minutes, ~1 GB RAM)
```

Reproduce the hardest-prime analysis (Section 5.4):

```bash
python verify_hardest_prime.py   # checks every claim for p = 3,807,728,761
```

Example output at 10^6:

```
Erdos-Straus Gateway Verification to 1,000,000
============================================================
Category                                            Count        %
-------------------------------------------------------------------
  p = 2                                                 1    0.001%
  p = 3 (mod 4)  [Prop. 4.1]                      39,322   50.093%
  p = 5 (mod 8)  [Prop. 4.2a]                     19,623   24.998%
  p = 17 (mod 24) [Prop. 4.2b]                     9,820   12.510%
  Case A  [Prop. 4.2c]                             5,192    6.614%
  Case B NQR7  [Prop. 4.3]                         2,271    2.893%
  Case B QR7  [Thm. 5.1]                           2,269    2.891%
-------------------------------------------------------------------
  PROVEN                                          78,498  100.0000%
  OPEN                                                 0    0.00000%

Max A needed: 79
ALL 78,498 PRIMES TO 1,000,000 VERIFIED.
```

For the full-scale segmented, checkpointed verification to 10^10 / 10^11 (16-worker pool, GPU sieve if CuPy is present):

```bash
python sessions/session20_corrected_10B.py 10000000000  results/session20_checkpoint.json   # 10^10
python sessions/session20_corrected_10B.py 100000000000 results/session21_checkpoint.json   # 10^11
```

Each run validates against three independent baselines (the full A-distribution at 10^6, the residual count at 10^7, and the residual count and max A at 10^9) before extending the range, and verifies every solution by direct evaluation of the identity. The `results/` directory holds the resulting checkpoints with per-decade milestone snapshots.

## Figures

Generate the paper figures (requires matplotlib):

```bash
python collect_stats.py     # regenerate stats_1M.json (10^6 sample)
python generate_figures.py  # fig1-fig5 as PDF and PNG
python generate_banner.py   # the README banner
```

| Figure | Description |
|--------|-------------|
| fig1_coverage | Prime classification hierarchy (pie charts) |
| fig2_A_distribution | Gateway parameter A distribution (bar + cumulative) |
| fig3_max_A | Max A growth across scales (10^3 to 10^11) |
| fig4_d_over_N | d/N ratio distribution showing the N^2 extension |
| fig5_gateway_diagram | Conceptual schematic of the decomposition |

## The Open Problem

Since every prime outside the residual class is covered algebraically, the conjecture reduces to:

> For every prime `p = 1 (mod 24)` that is a quadratic residue mod 7 with all prime factors of `(p+3)/4` congruent to `1 (mod 3)`, does there exist a bounded prime `A = 3 (mod 4)` such that `((p+A)/4)^2` has a divisor in the residue class `-p^2 * 4^{-1} (mod A)`?

This is a question about **equidistribution of divisors in residue classes** for structured integers, connected to work of Hooley and Tenenbaum.

## Data Provenance

The verification pipeline was corrected in July 2026: an earlier run of the 10^10 stage inverted the Case A / Case B classification and searched the wrong set of primes. All figures in this repository derive from the corrected, baseline-validated pipeline (`sessions/session20_corrected_10B.py`) and are cross-checked against the checkpoints in `results/`. Section 5 of the paper documents the correction.

## Citation

```bibtex
@article{erdos-straus-gateway-2026,
  title={Bounded Gateway Parameters for the {Erd\H{o}s--Straus} Conjecture},
  author={\'O Murch\'u, Macdara},
  year={2026},
  note={Preprint}
}
```

## License

MIT
