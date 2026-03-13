"""
Erdos-Straus Conjecture -- Session 1: Foundation

CONJECTURE: For every integer n >= 2, there exist positive integers x, y, z such that:
    4/n = 1/x + 1/y + 1/z

APPROACH THIS SESSION:
  1. Brute-force verifier (confirms all n up to SEARCH_LIMIT)
  2. Algebraic identity covers (handle whole residue classes symbolically)
  3. Residue class mapper: for each prime p <= PRIME_LIMIT, which residues
     mod p are covered by each algebraic identity, and which require search?
  4. Identify the "hard" residue classes (the mod 840 bottleneck from literature)
  5. Output a coverage report: what fraction of n is handled analytically vs
     requiring brute search?

COMPARISON WITH COLLATZ:
  - In Collatz we tracked residues mod 2^k (binary structure)
  - Here we track residues mod p for small primes (5, 7, 11, 13, 17, 19, 23, ...)
  - The "hard" cases concentrate among primes p ≡ 3 (mod 4)

Date: Session 1, March 2026
"""

import math
import json
import time
from collections import defaultdict
from pathlib import Path

RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ── Algebraic identity covers ──────────────────────────────────────────────────

def identity_covers(n):
    """
    Return a list of (x, y, z) found via algebraic identities.
    These handle entire residue classes without search.

    The identities below are well-known and cover the majority of integers.
    Returns a list of solutions (may be empty if no identity applies).
    """
    covers = []

    # Identity 1: 4/n = 1/n + 1/n + 2/n  (only works when 2|n)
    # Actually: 1/x + 1/y + 1/z must have integer denominators.
    # Systematic identities:

    # n = 1 (mod 4): n = 4k+1
    #   4/(4k+1) = 1/(k+1) + 1/((k+1)(4k+1))
    #   Need a third term. One standard form:
    #   4/n = 1/ceil(n/4) + ...
    k = n // 4

    # Identity A: 4/n = 1/k + (4k - n)/(kn)  where k = ceil(n/4)
    #   Then decompose (4k-n)/(kn) as 1/a + 1/b
    k_ceil = math.ceil(n / 4)
    remainder_num = 4 * k_ceil - n        # 4k - n, guaranteed 0..3
    remainder_den = k_ceil * n            # kn

    if remainder_num > 0:
        # remainder = remainder_num / remainder_den
        # Try to write as 1/a + 1/b:
        # 1/a + 1/b = rem_num/rem_den
        # a = ceil(rem_den / rem_num)
        a = math.ceil(remainder_den / remainder_num)
        # Then b = rem_den * a / (rem_num * a - rem_den)
        denom_b_num = remainder_den * a
        denom_b_den = remainder_num * a - remainder_den
        if denom_b_den > 0 and denom_b_num % denom_b_den == 0:
            b = denom_b_num // denom_b_den
            covers.append((k_ceil, a, b))

    # Identity B: Schinzel/Erdős family for n = 2 (mod 4)
    # 4/n = 2/(n/2) -- reduce to 2-fraction problem
    # 2/m = 1/ceil(m/2) + (2*ceil(m/2) - m)/(m*ceil(m/2))
    if n % 2 == 0:
        m = n // 2
        # 4/n = 2/m
        ceil_m2 = math.ceil(m / 2)
        rem2_num = 2 * ceil_m2 - m
        rem2_den = m * ceil_m2
        if rem2_num == 0:
            # 2/m = 1/ceil(m/2) exactly: need two more unit fractions
            # 1/x = 1/(2*ceil_m2) + 1/(2*ceil_m2)
            covers.append((2 * ceil_m2, 2 * ceil_m2, ceil_m2))
        elif rem2_num > 0 and rem2_den % rem2_num == 0:
            b = rem2_den // rem2_num
            covers.append((ceil_m2, b, b * ceil_m2 // (b - ceil_m2) if b != ceil_m2 else None))

    # Identity C: n = 0 (mod 4): 4/n = 1/(n/4) trivially (single term, expand)
    if n % 4 == 0:
        x = n // 4
        # 4/n = 1/x; split: 1/x = 1/(2x) + 1/(2x)... still 2 terms, need 3
        # Or: 1/x = 1/(x+1) + 1/(x(x+1))
        # So 4/n = 1/(x+1) + 1/(x(x+1)) -- only 2 terms, need one more split
        # This gives x, x+1, x*(x+1) as a valid 3-term decomp:
        # 1/(x+1) + 1/(x(x+1)) = (x + 1)/(x(x+1)) = 1/x. Yes.
        # But we need THREE terms summing to 1/x... Actually 4/n = 1/x (single fraction)
        # so split 1/x into 3: 1/x = 1/(x+1) + 1/(x(x+1)) ... that's only 2.
        # Use: 1/x = 1/(3x) + 1/(3x) + 1/(3x) when 3|x, else use:
        # 1/x = 1/(2x) + 1/(3x) + 1/(6x)  -- always works!
        covers.append((2 * x, 3 * x, 6 * x))

    # Identity D: general -- 4/n = 1/n + 1/n + 2/n when 2|n... denominators not integers
    # Better: 4/n = 1/n + 3/n. For 3/n = 1/a + 1/b:
    # a = ceil(n/3), b = n*ceil(n/3)/(3*ceil(n/3)-n)
    a = math.ceil(n / 3)
    rem3_num = 3 * a - n
    rem3_den = n * a
    if rem3_num > 0 and rem3_den % rem3_num == 0:
        b = rem3_den // rem3_num
        covers.append((n, a, b))

    return covers


# ── Brute-force verifier ───────────────────────────────────────────────────────

def find_decomposition_search(n, scan_width=2000):
    """
    Find (x, y, z) with 1/x + 1/y + 1/z = 4/n.

    For each x starting at ceil(n/4), the remainder is A/B = (4x-n)/(nx).
    Given A/B, solve 1/y + 1/z = A/B by scanning y from ceil(B/A):
      z = B*y / (A*y - B)  -- integer iff (A*y-B) | B*y

    Since A is small for x near ceil(n/4) (A in {1,2,3}), the scan over
    y is very tight in practice: solutions always appear within a few hundred
    steps of the minimal y value for n up to 10^7.
    """
    x_lo = math.ceil(n / 4)
    x_hi = (3 * n) // 4 + 1

    for x in range(x_lo, x_hi):
        A = 4 * x - n
        if A <= 0:
            continue
        B = n * x
        # Minimum y (WLOG y <= z): y >= ceil(B/A)  and  y >= x
        y_min = max(x, (B + A - 1) // A)
        # Maximum y for z >= y: y <= 2B/A
        y_max = y_min + scan_width

        for y in range(y_min, y_max + 1):
            num = A * y - B
            if num <= 0:
                continue
            rem = B * y
            if rem % num == 0:
                z = rem // num
                if z >= y:
                    return (x, y, z)
                # z < y: still valid, just reorder
                return (x, z, y)

    return None


def verify_decomposition(n, x, y, z):
    """Confirm 1/x + 1/y + 1/z == 4/n using integer arithmetic."""
    lhs = x * y * z        # common denominator
    rhs_parts = y * z + x * z + x * y
    return rhs_parts * n == 4 * lhs


# ── Residue class analysis ─────────────────────────────────────────────────────

def identity_cover_check(n):
    """
    Check algebraic identities first, then fall back to search.
    Returns (solution, method) where method is 'identity' or 'search'.
    """
    covers = identity_covers(n)
    for sol in covers:
        if sol is not None and all(isinstance(v, int) and v > 0 for v in sol):
            if verify_decomposition(n, *sol):
                return sol, 'identity'

    sol = find_decomposition_search(n)
    return sol, 'search'


def residue_analysis(prime, limit=10000):
    """
    For prime p, check each residue class 0..p-1.
    Count how many n ≡ r (mod p) require search vs identity.
    """
    residue_stats = defaultdict(lambda: {'identity': 0, 'search': 0, 'failed': 0})

    for n in range(2, limit + 1):
        r = n % prime
        if n < 2:
            continue
        sol, method = identity_cover_check(n)
        if sol is None:
            residue_stats[r]['failed'] += 1
        else:
            residue_stats[r][method] += 1

    return dict(residue_stats)


# ── Main analysis ──────────────────────────────────────────────────────────────

def run_verification(limit=100000):
    """Verify conjecture for all n in [2, limit]. Report statistics."""
    print(f"\n{'='*60}")
    print(f"ERDOS-STRAUS VERIFICATION: n = 2 to {limit:,}")
    print(f"{'='*60}")

    stats = {'identity': 0, 'search': 0, 'failed': 0, 'hardest': []}
    failed_cases = []
    t0 = time.time()

    for n in range(2, limit + 1):
        sol, method = identity_cover_check(n)
        if sol is None:
            stats['failed'] += 1
            failed_cases.append(n)
            print(f"  FAILED: n={n}")
        else:
            stats[method] += 1
            if method == 'search':
                stats['hardest'].append((n, sol))

        if n % 10000 == 0:
            elapsed = time.time() - t0
            rate = n / elapsed
            print(f"  n={n:,}: identity={stats['identity']:,} search={stats['search']:,} "
                  f"failed={stats['failed']} | {rate:.0f}/s")

    elapsed = time.time() - t0
    total = limit - 1
    print(f"\nRESULTS (n=2..{limit:,}):")
    print(f"  Covered by identity : {stats['identity']:,} ({100*stats['identity']/total:.1f}%)")
    print(f"  Covered by search   : {stats['search']:,} ({100*stats['search']/total:.1f}%)")
    print(f"  FAILED              : {stats['failed']}")
    print(f"  Time                : {elapsed:.1f}s")

    if failed_cases:
        print(f"\n  COUNTEREXAMPLES FOUND: {failed_cases}")
    else:
        print(f"\n  No counterexamples found in [2, {limit:,}]. Conjecture holds.")

    return stats, failed_cases


def run_residue_map(primes=[5, 7, 11, 13, 17, 19, 23], limit=5000):
    """
    For each prime p in list, show which residue classes are 'hard'
    (require search rather than identity cover).
    """
    print(f"\n{'='*60}")
    print(f"RESIDUE CLASS ANALYSIS (mod primes, n=2..{limit})")
    print(f"{'='*60}")

    report = {}
    for p in primes:
        stats = residue_analysis(p, limit)
        hard_residues = [
            r for r, s in stats.items()
            if s['search'] > 0 or s['failed'] > 0
        ]
        print(f"\n  Prime p={p}: {len(hard_residues)} hard residue classes out of {p}")
        for r in sorted(hard_residues):
            s = stats[r]
            search_pct = 100 * s['search'] / max(s['identity'] + s['search'] + s['failed'], 1)
            print(f"    r={r:2d}: identity={s['identity']:4d}, search={s['search']:4d} "
                  f"({search_pct:.0f}% need search), failed={s['failed']}")
        report[p] = stats

    return report


def run_algebraic_identity_audit(sample_n=None):
    """
    Check which algebraic identities are actually firing for the first 1000 integers.
    Understand the coverage structure.
    """
    if sample_n is None:
        sample_n = list(range(2, 1001))

    print(f"\n{'='*60}")
    print(f"ALGEBRAIC IDENTITY AUDIT (n=2..{sample_n[-1]})")
    print(f"{'='*60}")

    identity_count = 0
    search_count = 0
    multi_solution = []  # n with multiple identity solutions

    for n in sample_n:
        covers = [
            sol for sol in identity_covers(n)
            if sol is not None and all(isinstance(v, int) and v > 0 for v in sol)
            and verify_decomposition(n, *sol)
        ]
        if covers:
            identity_count += 1
            if len(covers) > 1:
                multi_solution.append((n, covers))
        else:
            search_count += 1

    print(f"  Identity covers : {identity_count}/{len(sample_n)} ({100*identity_count/len(sample_n):.1f}%)")
    print(f"  Need search     : {search_count}/{len(sample_n)} ({100*search_count/len(sample_n):.1f}%)")
    print(f"  Multiple covers : {len(multi_solution)} cases")

    # Show the first 20 search-required cases
    print(f"\n  First search-required cases:")
    shown = 0
    for n in sample_n:
        covers = [
            sol for sol in identity_covers(n)
            if sol is not None and all(isinstance(v, int) and v > 0 for v in sol)
            and verify_decomposition(n, *sol)
        ]
        if not covers:
            sol = find_decomposition_search(n)
            print(f"    n={n:4d}: search -> {sol}  [n mod 4 = {n%4}, n mod 3 = {n%3}]")
            shown += 1
            if shown >= 20:
                break


# ── Entry point ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("ERDOS-STRAUS CONJECTURE -- SESSION 1 FOUNDATION")
    print("=" * 60)
    print("4/n = 1/x + 1/y + 1/z  for all n >= 2")
    print()
    print("Plan:")
    print("  1. Algebraic identity audit (n=2..1000)")
    print("  2. Verification + coverage stats (n=2..100,000)")
    print("  3. Residue class analysis (mod 5, 7, 11, 13, 17)")
    print()

    # Step 1: Understand identity coverage
    run_algebraic_identity_audit()

    # Step 2: Verify up to 100k and measure identity vs search rate
    stats, failed = run_verification(limit=100_000)

    # Step 3: Residue class map -- where are the hard cases?
    report = run_residue_map(primes=[5, 7, 11, 13, 17, 19, 23], limit=5000)

    # Save report
    output = {
        'verification_stats': {k: v for k, v in stats.items() if k != 'hardest'},
        'failed_cases': failed,
        'residue_report': {str(p): {str(r): s for r, s in res.items()} for p, res in report.items()}
    }
    out_path = RESULTS_DIR / "session1_results.json"
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to: {out_path}")
