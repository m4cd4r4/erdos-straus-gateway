# A Finite Algebraic Covering System for the Erdos-Straus Conjecture to 10^9

**Gateway Decompositions via Divisors of N^2**

## Overview

The [Erdos-Straus conjecture](https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93Straus_conjecture) (1948) asserts that for every integer n >= 2, the equation

    4/n = 1/x + 1/y + 1/z

has a solution in positive integers x, y, z.

This repository contains the paper, verification code, and figures for a unified algebraic approach that covers **all 50,847,534 primes up to 10^9** using just **24 values** of a single auxiliary parameter A (all prime, all <= 239).

## Key Result

**Theorem.** For every prime p <= 10^9, there exists a prime A <= 239 with A = 3 (mod 4) and a divisor d of ((p+A)/4)^2 satisfying d = -p^2 * 4^{-1} (mod A), such that

    x = (p+A)/4,  y = (pN+d)/A,  z = pNy/d

are positive integers giving 4/p = 1/x + 1/y + 1/z.

The critical insight is that the integrality condition requires d | N^2 (where N = (p+A)/4), which is **strictly weaker** than d | N. About 1.8% of the hardest primes require this extension.

## Paper

The paper (`erdos_straus_gateway.tex`) contains:

- **Section 3**: The Gateway Decomposition theorem and the d | N^2 lemma
- **Section 4**: Algebraic existence proofs covering ~97.6% of all primes
- **Section 5**: Computational verification for the remaining ~2.4% (Case B QR7 primes)
- **Section 6**: Discussion of the bounded-A phenomenon and the path to a full proof

Compiled PDF: `erdos_straus_gateway.pdf`

## Verification

Run the self-contained verification script:

```bash
# Verify all primes to 10^6 (~10 seconds)
python verify.py

# Verify to 10^7 (~2 minutes)
python verify.py 10000000

# Verify to 10^8 (~30 minutes)
python verify.py 100000000

# Verify to 10^9 (~90 minutes, ~1GB RAM)
python verify.py 1000000000
```

No dependencies beyond Python 3.6+ standard library.

Example output at 10^6:

```
Erdos-Straus Gateway Verification to 1,000,000
============================================================
Category                                            Count        %
-------------------------------------------------------------------
  p = 2                                                 1    0.001%
  p = 3 (mod 4)  [Prop. 4.1]                      39,322   50.093%
  p = 5 (mod 8)  [Prop. 4.2a]                     19,617   24.990%
  p = 17 (mod 24) [Prop. 4.2b]                     9,828   12.521%
  Case A  [Prop. 4.2c]                              5,192    6.614%
  Case B NQR7  [Prop. 4.3]                          2,271    2.893%
  Case B QR7  [Thm. 5.1]                            2,267    2.889%
-------------------------------------------------------------------
  PROVEN                                           78,498  100.0000%
  OPEN                                                  0    0.00000%

ALL 78,498 PRIMES TO 1,000,000 VERIFIED.
```

## Figures

Generate the paper figures (requires matplotlib):

```bash
python collect_stats.py    # Generates stats_1M.json
python generate_figures.py # Generates fig1-fig5 as PDF and PNG
```

| Figure | Description |
|--------|-------------|
| fig1_coverage | Prime classification hierarchy (pie charts) |
| fig2_A_distribution | Gateway parameter A distribution (bar + cumulative) |
| fig3_max_A | Max A stabilisation across scales (10^3 to 10^9) |
| fig4_d_over_N | d/N ratio distribution showing N^2 extension |
| fig5_gateway_diagram | Conceptual schematic of the decomposition |

## The Open Problem

The conjecture reduces to:

> For every prime p = 1 (mod 24) that is a quadratic residue mod 7 with all prime factors of (p+3)/4 congruent to 1 (mod 3), does there exist a bounded prime A = 3 (mod 4) such that ((p+A)/4)^2 has a divisor in the residue class -p^2 * 4^{-1} (mod A)?

This is a question about **equidistribution of divisors in residue classes** for structured integers, connected to work of Hooley and Tenenbaum.

## Citation

If you use this work, please cite:

```bibtex
@article{erdos-straus-gateway-2026,
  title={A Finite Algebraic Covering System for the {Erd\H{o}s--Straus} Conjecture to $10^9$},
  author={Macdara},
  year={2026},
  note={Preprint}
}
```

## License

MIT
