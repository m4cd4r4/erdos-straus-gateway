"""Calibration: cost of Tier 2 - independently reclassifying EVERY prime in
[10^11, 1.212x10^12], not just the A>=100 witnesses Tier 1 covered.

This is a different regime from Tier 1. Tier 1 touched 19,540 primes (the
expensive tail). Tier 2 means classifying all ~33.5 BILLION primes in
[10^11, 10^12] alone (real count, from the checkpoint's own milestone data -
not an estimate), matching the <=10^11 table's actual standard: every prime,
not a sample.

Independent implementation (sympy, not the primary engine's hand-rolled
trial-division + Pollard-Brent) of the FULL cascade:
  p = 2                                   -> p2
  p = 3 mod 4                             -> thm1
  p = 5 mod 8                             -> thm2
  p = 17 mod 24                           -> thm4
  p = 1 mod 24, some factor of (p+3)/4
      not = 1 mod 3                       -> caseA
  p = 1 mod 24, Case B, NQR mod 7         -> nqr7
  otherwise (Case B, QR mod 7)            -> gateway search (the expensive one)

Two calibration windows at different scales (just above 1e11, just below
1e12) since factorisation cost - and therefore total cost - may not be flat
across the range. Extrapolates using the REAL prime-count difference between
milestones (33,489,857,205 in [1e11,1e12]), not a density formula.

Run: python calibrate_tier2_verify.py [WIDTH]
"""
import sys, time, json
from sympy import factorint, primerange, isprime, legendre_symbol

WIDTH = int(sys.argv[1]) if len(sys.argv) > 1 else 2_000_000
CKPT = 'results/session23_checkpoint.json'

WINDOWS = [
    ('just above 1e11', 100_000_000_000 + 1),
    ('just below 1e12', 1_000_000_000_000 - WIDTH),
]

GW = [q for q in primerange(3, 20000) if q % 4 == 3]


def divisor_residues(facts, B):
    D = {1 % B}
    for q, e in facts.items():
        if q % B == 0:
            continue
        qm, pw, new = q % B, 1, set()
        for j in range(2 * e + 1):
            if j:
                pw = (pw * qm) % B
            for r in D:
                new.add((r * pw) % B)
        D = new
    return D


def classify(p):
    """Returns a stage tag; does the full gateway search only when needed."""
    if p % 4 == 3:
        return 'thm1'
    if p % 8 == 5:
        return 'thm2'
    if p % 24 == 17:
        return 'thm4'
    # p = 1 mod 24
    m = (p + 3) // 4
    if any(q % 3 != 1 for q in factorint(m)):
        return 'caseA'
    if legendre_symbol(p, 7) == -1:
        return 'nqr7'
    # Case B, QR mod 7 - the expensive gateway search
    for A in GW:
        N = (p + A) // 4
        if N % A == 0:
            continue
        t = (-4 * N * N) % A
        if t in divisor_residues(factorint(N), A):
            return f'gw:{A}'
    return 'OPEN'


def main():
    d = json.load(open(CKPT))
    p_1e11 = d['milestones']['100000000000']['primes']
    p_1e12 = d['milestones']['1000000000000']['primes']
    total_primes_1e11_1e12 = p_1e12 - p_1e11
    print(f'REAL prime count to classify in [1e11, 1e12]: {total_primes_1e11_1e12:,}')
    print(f'(window width for this calibration: {WIDTH:,})')
    print()

    rates = []
    for label, lo in WINDOWS:
        hi = lo + WIDTH
        primes = list(primerange(lo, hi))
        print(f'{label}: [{lo:,}, {hi:,}), {len(primes):,} primes')

        stage_counts = {}
        t0 = time.time()
        for p in primes:
            s = classify(p)
            key = 'gw' if s.startswith('gw:') else s
            stage_counts[key] = stage_counts.get(key, 0) + 1
        elapsed = time.time() - t0

        rate = len(primes) / elapsed if elapsed > 0 else float('inf')
        rates.append((label, lo, len(primes), elapsed, rate))
        print(f'  {elapsed:.2f}s, {rate:.0f} primes/s')
        print(f'  stages: {stage_counts}')
        print()

    print('=' * 70)
    print('EXTRAPOLATION to the full [1e11, 1e12] range '
          f'({total_primes_1e11_1e12:,} real primes):')
    for label, lo, n, elapsed, rate in rates:
        proj_s = total_primes_1e11_1e12 / rate
        print(f'  using the "{label}" rate ({rate:.0f} primes/s):')
        print(f'    {proj_s:,.0f} s = {proj_s/3600:,.1f} core-hours '
              f'= {proj_s/3600/24:,.1f} core-days')
        for workers in (10, 16):
            print(f'      wall time at {workers} workers: '
                  f'{proj_s/3600/workers/24:,.2f} days')
    print()
    print('NOTE: this is [1e11,1e12] only. The actual endpoint is 1.212e12;')
    print('scale up by roughly another 21% for the full extended-search range.')

    json.dump(dict(width=WIDTH, windows=[(l, lo, n, e, r) for l, lo, n, e, r in rates],
                    total_primes_1e11_1e12=total_primes_1e11_1e12),
              open('results/calibrate_tier2_results.json', 'w'))


if __name__ == '__main__':
    main()
