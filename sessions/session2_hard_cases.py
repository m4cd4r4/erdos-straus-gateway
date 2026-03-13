"""
Erdos-Straus Conjecture -- Session 2: Hard Case Analysis

SESSION 1 FINDING: Residue classes mod single primes (5,7,11,...,23) show
UNIFORM ~37-38% search rate across ALL residues. No single prime creates a
special hard class. The hard structure must be in COMBINED residues.

THIS SESSION:
  1. Prime-only analysis: for primes p, measure the "x-distance" from ceil(p/4)
     to the first solution. Distribution by p mod 4 and p mod 8.
  2. CRT analysis: mod 60 and mod 840 combined residues.
     Find which (n mod 840) values have NO algebraic identity cover.
  3. Search-distance distribution: histogram of how far x must go.
  4. Identify algebraic structure of hard primes p = 3 (mod 4).

GOAL: Quantify the gap between "identity covers" and "all n".
      Find if any combined residue class has systematically harder solutions.

Date: Session 2, March 2026
"""

import math
import json
import time
from collections import defaultdict
from pathlib import Path

RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)


# ── Sieve and core utilities ───────────────────────────────────────────────────

def sieve_primes(limit):
    """Sieve of Eratosthenes."""
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, limit + 1) if is_prime[i]]


def find_solution_with_distance(n, scan_width=5000):
    """
    Find (x, y, z) for 4/n = 1/x + 1/y + 1/z.
    Returns (x, y, z, x_dist, y_dist) where:
      x_dist = x - ceil(n/4)   (steps x is above minimal)
      y_dist = y - ceil(B/A)   (steps y is above minimal for that x)
    """
    x_lo = math.ceil(n / 4)
    x_hi = (3 * n) // 4 + 1

    for x in range(x_lo, min(x_lo + scan_width, x_hi)):
        A = 4 * x - n
        if A <= 0:
            continue
        B = n * x
        y_min = max(x, (B + A - 1) // A)

        for y in range(y_min, y_min + scan_width + 1):
            num = A * y - B
            if num <= 0:
                continue
            rem = B * y
            if rem % num == 0:
                z = rem // num
                if z > 0:
                    return (x, y, z, x - x_lo, y - y_min)

    return None


def verify(n, x, y, z):
    return (y*z + x*z + x*y) * n == 4 * x * y * z


# ── Analysis 1: Prime-only x-distance distribution ────────────────────────────

def prime_distance_analysis(prime_limit=50000):
    """
    For each prime p up to prime_limit, find the solution and record:
    - x_dist: how far x is from ceil(p/4)
    - p mod 4 (1 or 3 -- primes > 2 are always odd)
    - p mod 8 (1, 3, 5, 7)
    """
    print(f"\n{'='*60}")
    print(f"PRIME-ONLY X-DISTANCE ANALYSIS (primes up to {prime_limit:,})")
    print(f"{'='*60}")

    primes = sieve_primes(prime_limit)
    t0 = time.time()

    dist_by_mod4 = defaultdict(list)   # mod4 -> [x_dist, ...]
    dist_by_mod8 = defaultdict(list)
    hardest = []   # (x_dist, p, solution)
    failed = []

    for p in primes:
        result = find_solution_with_distance(p)
        if result is None:
            failed.append(p)
            continue
        x, y, z, x_dist, y_dist = result
        dist_by_mod4[p % 4].append(x_dist)
        dist_by_mod8[p % 8].append(x_dist)
        if x_dist >= 5:
            hardest.append((x_dist, p, (x, y, z)))

    elapsed = time.time() - t0
    print(f"  Primes checked: {len(primes):,} in {elapsed:.1f}s")
    print(f"  Failed: {failed}")

    print(f"\n  X-DISTANCE BY p mod 4:")
    for mod, dists in sorted(dist_by_mod4.items()):
        avg = sum(dists) / len(dists)
        max_d = max(dists)
        pct_zero = 100 * sum(1 for d in dists if d == 0) / len(dists)
        print(f"    p = {mod} (mod 4): n={len(dists):,}, avg={avg:.3f}, "
              f"max={max_d}, x_dist=0: {pct_zero:.1f}%")

    print(f"\n  X-DISTANCE BY p mod 8:")
    for mod, dists in sorted(dist_by_mod8.items()):
        avg = sum(dists) / len(dists)
        max_d = max(dists)
        pct_zero = 100 * sum(1 for d in dists if d == 0) / len(dists)
        print(f"    p = {mod} (mod 8): n={len(dists):,}, avg={avg:.3f}, "
              f"max={max_d}, x_dist=0: {pct_zero:.1f}%")

    print(f"\n  HARDEST PRIMES (x_dist >= 5):")
    hardest.sort(reverse=True)
    for x_dist, p, sol in hardest[:20]:
        print(f"    p={p:7d} (mod 4={p%4}, mod 8={p%8}): x_dist={x_dist}, sol={sol}")

    # Build histogram
    all_dists = [d for dlist in dist_by_mod4.values() for d in dlist]
    hist = defaultdict(int)
    for d in all_dists:
        hist[min(d, 20)] += 1

    print(f"\n  X-DISTANCE HISTOGRAM (all primes):")
    total = len(all_dists)
    cumulative = 0
    for d in range(21):
        count = hist[d]
        cumulative += count
        label = f"{d}" if d < 20 else ">=20"
        bar = '#' * (count * 40 // max(hist.values()))
        print(f"    dist={label:4s}: {count:6,} ({100*count/total:5.1f}%) "
              f"cumul={100*cumulative/total:5.1f}%  {bar}")

    return {
        'by_mod4': {k: {'avg': sum(v)/len(v), 'max': max(v), 'n': len(v)} for k,v in dist_by_mod4.items()},
        'by_mod8': {k: {'avg': sum(v)/len(v), 'max': max(v), 'n': len(v)} for k,v in dist_by_mod8.items()},
        'hardest': [(d, p) for d, p, _ in hardest[:50]],
        'failed': failed
    }


# ── Analysis 2: CRT combined residue analysis ──────────────────────────────────

def crt_residue_analysis(modulus=840, limit=100000):
    """
    For each combined residue class r mod 'modulus', measure:
    - fraction covered by algebraic identities
    - average x_dist for those needing search
    - max x_dist

    mod 840 = mod (8 * 3 * 5 * 7) covers the main algebraic identity structure.
    The literature identifies hard cases concentrate at specific mod-840 residues.
    """
    print(f"\n{'='*60}")
    print(f"CRT RESIDUE ANALYSIS (mod {modulus}, n=2..{limit:,})")
    print(f"{'='*60}")

    # Quick algebraic identity check (from session 1)
    def has_identity(n):
        k = n // 4
        k_ceil = math.ceil(n / 4)
        rem_num = 4 * k_ceil - n
        rem_den = k_ceil * n
        if rem_num > 0:
            a = math.ceil(rem_den / rem_num)
            denom_b_num = rem_den * a
            denom_b_den = rem_num * a - rem_den
            if denom_b_den > 0 and denom_b_num % denom_b_den == 0:
                return True
        if n % 4 == 0:
            return True  # identity C always works
        a = math.ceil(n / 3)
        rem3_num = 3 * a - n
        rem3_den = n * a
        if rem3_num > 0 and rem3_den % rem3_num == 0:
            return True
        return False

    residue_stats = defaultdict(lambda: {
        'identity': 0, 'search': 0, 'x_dist_sum': 0,
        'x_dist_max': 0, 'y_dist_sum': 0
    })

    t0 = time.time()
    for n in range(2, limit + 1):
        r = n % modulus
        if has_identity(n):
            residue_stats[r]['identity'] += 1
        else:
            result = find_solution_with_distance(n)
            if result:
                _, _, _, x_dist, y_dist = result
                residue_stats[r]['search'] += 1
                residue_stats[r]['x_dist_sum'] += x_dist
                residue_stats[r]['x_dist_max'] = max(residue_stats[r]['x_dist_max'], x_dist)
                residue_stats[r]['y_dist_sum'] += y_dist

    elapsed = time.time() - t0
    print(f"  Processed {limit:,} values in {elapsed:.1f}s")

    # Find the hardest residue classes (highest avg x_dist)
    hard_residues = []
    for r, s in residue_stats.items():
        total = s['identity'] + s['search']
        if s['search'] > 0:
            avg_xd = s['x_dist_sum'] / s['search']
            identity_pct = 100 * s['identity'] / total
            hard_residues.append((avg_xd, r, s['x_dist_max'], identity_pct, s['search']))

    hard_residues.sort(reverse=True)

    print(f"\n  TOP 30 HARDEST RESIDUE CLASSES (mod {modulus}) by avg x-distance:")
    print(f"  {'r':>6}  {'r mod 4':>7}  {'r mod 3':>7}  {'avg_xd':>7}  "
          f"{'max_xd':>7}  {'id_pct':>7}  {'n_search':>9}")
    for avg_xd, r, max_xd, id_pct, n_search in hard_residues[:30]:
        print(f"  {r:6d}  {r%4:7d}  {r%3:7d}  {avg_xd:7.3f}  "
              f"{max_xd:7d}  {id_pct:6.1f}%  {n_search:9,}")

    # Which mod-840 classes have ZERO identity coverage?
    zero_id = [(r, s) for r, s in residue_stats.items() if s['identity'] == 0]
    print(f"\n  Residue classes with ZERO identity coverage: {len(zero_id)}")
    for r, s in sorted(zero_id)[:20]:
        print(f"    r={r:4d} (mod 4={r%4}, mod 3={r%3}, mod 5={r%5}, mod 7={r%7}): "
              f"all {s['search']} need search, max_xd={s['x_dist_max']}")

    return {
        'top_hard': [(avg_xd, r, max_xd, id_pct) for avg_xd, r, max_xd, id_pct, _ in hard_residues[:50]],
        'zero_identity_residues': [r for r, _ in zero_id]
    }


# ── Analysis 3: Prime gap -- hard primes pattern ───────────────────────────────

def hard_prime_pattern(prime_limit=200000):
    """
    Find the hardest primes (largest x_dist) and look for patterns:
    - Are they clustered in specific mod-24 or mod-120 residue classes?
    - Do they have specific small prime factor structure?
    - Are they related to Mersenne or other special prime families?
    """
    print(f"\n{'='*60}")
    print(f"HARD PRIME PATTERN ANALYSIS (primes up to {prime_limit:,})")
    print(f"{'='*60}")

    primes = sieve_primes(prime_limit)
    t0 = time.time()

    results = []
    for p in primes:
        result = find_solution_with_distance(p, scan_width=200)
        if result:
            x, y, z, x_dist, y_dist = result
            results.append((x_dist, p, x, y, z))

    results.sort(reverse=True)
    elapsed = time.time() - t0
    print(f"  {len(primes):,} primes checked in {elapsed:.1f}s")

    print(f"\n  TOP 40 HARDEST PRIMES:")
    print(f"  {'p':>10}  {'mod4':>5}  {'mod8':>5}  {'mod12':>6}  "
          f"{'mod24':>6}  {'xd':>5}  {'x':>12}  sol_structure")
    for x_dist, p, x, y, z in results[:40]:
        x_lo = math.ceil(p / 4)
        print(f"  {p:10d}  {p%4:5d}  {p%8:5d}  {p%12:6d}  "
              f"{p%24:6d}  {x_dist:5d}  {x:12d}  "
              f"x={x_lo}+{x_dist}, y_ratio={y//x if x else 0}")

    # Histogram of x_dist for primes (1 only)
    p3mod4 = [r for r in results if r[1] % 4 == 3]
    p1mod4 = [r for r in results if r[1] % 4 == 1]

    print(f"\n  X-DIST COMPARISON:")
    print(f"  p=3 (mod 4): n={len(p3mod4):,}, "
          f"avg={sum(r[0] for r in p3mod4)/max(len(p3mod4),1):.4f}, "
          f"max={max((r[0] for r in p3mod4), default=0)}")
    print(f"  p=1 (mod 4): n={len(p1mod4):,}, "
          f"avg={sum(r[0] for r in p1mod4)/max(len(p1mod4),1):.4f}, "
          f"max={max((r[0] for r in p1mod4), default=0)}")

    # Distribution of x_dist=0 (solution at minimal x)
    xd0_3mod4 = sum(1 for r in p3mod4 if r[0] == 0)
    xd0_1mod4 = sum(1 for r in p1mod4 if r[0] == 0)
    print(f"\n  x_dist=0 rate (solution at minimal x=ceil(p/4)):")
    print(f"  p=3 (mod 4): {100*xd0_3mod4/max(len(p3mod4),1):.1f}%")
    print(f"  p=1 (mod 4): {100*xd0_1mod4/max(len(p1mod4),1):.1f}%")

    return [(x_dist, p) for x_dist, p, *_ in results[:100]]


# ── Entry point ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("ERDOS-STRAUS -- SESSION 2: HARD CASE ANALYSIS")
    print("=" * 60)
    print()
    print("Session 1 established: conjecture holds for n=2..100,000.")
    print("Session 1 finding: no single prime creates a structurally hard residue class.")
    print("Session 2 goal: find WHERE the hard structure actually lives.")
    print()

    all_results = {}

    # Analysis 1: Prime x-distance distribution
    r1 = prime_distance_analysis(prime_limit=50000)
    all_results['prime_distance'] = r1

    # Analysis 2: CRT mod-840 hard residue classes
    r2 = crt_residue_analysis(modulus=840, limit=50000)
    all_results['crt_840'] = r2

    # Analysis 3: Hard prime patterns
    r3 = hard_prime_pattern(prime_limit=100000)
    all_results['hard_prime_pattern'] = [(d, p) for d, p in r3]

    # Save
    out_path = RESULTS_DIR / "session2_results.json"
    with open(out_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to: {out_path}")
    print("\nSESSION 2 COMPLETE.")
    print("Next: session3_algebraic_covers.py -- extend identity coverage to eliminate hard classes")
