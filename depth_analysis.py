#!/usr/bin/env python3
"""
Depth analysis for the gateway decomposition.

For each Case B QR7 prime p and its successful gateway (A, d):
  - Find the NQR prime factor q of N_A = (p+A)/4 that "saves" the gateway
  - Record v_q(N_A) = NQR depth
  - Check whether the A=7 sharp theorem holds (depth 1 always sufficient for A=7)
  - Check whether depth >= 2 is ever NEEDED (i.e., depth=1 saves us or not)

Key question (Q4 from companion paper):
  Does NQR depth >= 2 UNCONDITIONALLY guarantee gateway success?
  i.e., if v_q(N_A) >= 2 for some NQR prime q, does the gateway always succeed?

Usage:
  python depth_analysis.py           # analyse to 10^6
  python depth_analysis.py 10000000  # analyse to 10^7
"""

import sys
from math import isqrt

LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000

# === Utilities ===

def sieve_primes(n):
    ip = bytearray([1]) * (n + 1)
    ip[0] = ip[1] = 0
    for i in range(2, isqrt(n) + 1):
        if ip[i]:
            ip[i*i::i] = bytearray(len(ip[i*i::i]))
    return ip

def factorise(n):
    """Return dict {p: e} for n = prod p^e."""
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors

def legendre(a, p):
    """Legendre symbol (a/p) for odd prime p."""
    a = a % p
    if a == 0:
        return 0
    r = pow(a, (p - 1) // 2, p)
    return -1 if r == p - 1 else r

def divisors_of_nsq(n):
    """Return all divisors of n^2."""
    factors = factorise(n)
    divs = [1]
    for p, e in factors.items():
        new_divs = []
        pe = 1
        for k in range(2 * e + 1):
            new_divs.extend(d * pe for d in divs)
            pe *= p
        divs = new_divs
    return divs

def gateway_works(p, A):
    """
    Try gateway with parameter A for prime p.
    Returns (d, N) if it works, None otherwise.
    Requires A odd, A = 3 mod 4, 4 | (p+A).
    """
    if (p + A) % 4 != 0:
        return None
    N = (p + A) // 4
    B = p * N
    target = (-B) % A  # need d ≡ -B ≡ -4N² (mod A)
    for d in divisors_of_nsq(N):
        if d % A == target:
            # Verify integrality of z
            if (B * ((B + d) // A)) % d == 0:
                return (d, N)
    return None

def nqr_depth(N, A):
    """
    Return (min_nqr_valuation, has_nqr_factor):
      - has_nqr_factor: True if N has any prime factor q with (q/A) = -1
      - min_nqr_valuation: minimum v_q(N) among NQR prime factors q
    """
    factors = factorise(N)
    nqr_vals = []
    for q, e in factors.items():
        if legendre(q, A) == -1:
            nqr_vals.append(e)
    if not nqr_vals:
        return (None, False)
    return (min(nqr_vals), True)

def is_case_b_qr7(p):
    """Return True if p is a Case B QR7 prime."""
    if p % 24 != 1:
        return False
    # p must be QR mod 7: p ≡ 1, 2, or 4 mod 7
    if p % 7 not in (1, 2, 4):
        return False
    # All prime factors of (p+3)/4 must be ≡ 1 mod 3
    n = (p + 3) // 4
    factors = factorise(n)
    for q in factors:
        if q % 3 != 1:
            return False
    return True

# Candidate A values (primes ≡ 3 mod 4, up to 239)
CANDIDATE_AS = [
    a for a in range(3, 240)
    if a % 4 == 3 and all(a % d != 0 for d in range(2, isqrt(a)+1))
]

# === Main analysis ===

print(f"Depth Analysis for Case B QR7 primes up to {LIMIT:,}")
print("=" * 65)

ip = sieve_primes(LIMIT)
primes = [p for p in range(5, LIMIT+1) if ip[p]]

# Statistics
total_case_b = 0
solved_depth1 = 0      # saved by a depth-1 NQR factor
solved_depth2plus = 0  # saved by a depth->=2 NQR factor
solved_no_nqr_needed = 0  # this shouldn't happen (A fails means QA-smooth)

# For A=7 specifically: is depth-1 always sufficient?
a7_failures_with_nqr = 0  # A=7 failed despite having NQR factor

# Hard cases: primes where all depth-1 attempts fail and depth-2 saves
depth2_saves = []

# Per-A analysis
a_depth_stats = {}  # A -> {depth: count}

for p in primes:
    if not is_case_b_qr7(p):
        continue
    total_case_b += 1

    # Find the successful A and analyse depth
    saved = False
    for A in CANDIDATE_AS:
        if (p + A) % 4 != 0:
            continue
        result = gateway_works(p, A)
        if result is not None:
            d, N = result
            depth, has_nqr = nqr_depth(N, A)
            if not has_nqr:
                # Shouldn't happen: Q_A-smooth N should fail
                solved_no_nqr_needed += 1
            elif depth == 1:
                solved_depth1 += 1
            else:
                solved_depth2plus += 1
                depth2_saves.append((p, A, depth))

            # Per-A stats
            if A not in a_depth_stats:
                a_depth_stats[A] = {}
            a_depth_stats[A][depth if has_nqr else 0] = \
                a_depth_stats[A].get(depth if has_nqr else 0, 0) + 1

            saved = True
            break

    if not saved:
        print(f"  WARNING: p={p} not solved!")

# === A=7 sharp theorem check ===
print("\n--- A=7 Sharp Theorem Verification ---")
print("(Gateway A=7 fails iff N_7 is Q_7-smooth,")
print(" succeeds iff N_7 has ANY NQR factor mod 7)")
a7_with_nqr_fails = 0
a7_total = 0
for p in primes:
    if not is_case_b_qr7(p):
        continue
    if (p + 7) % 4 != 0:
        continue
    a7_total += 1
    N7 = (p + 7) // 4
    _, has_nqr = nqr_depth(N7, 7)
    result = gateway_works(p, 7)
    if has_nqr and result is None:
        a7_with_nqr_fails += 1
        print(f"  p={p}: N_7={N7} has NQR factor but A=7 FAILED (unexpected!)")
    if not has_nqr and result is not None:
        print(f"  p={p}: N_7={N7} is Q_7-smooth but A=7 SUCCEEDED (impossible!)")

print(f"  Checked {a7_total} Case B QR7 primes with 4|(p+7)")
if a7_with_nqr_fails == 0:
    print("  CONFIRMED: A=7 sharp theorem holds -")
    print("  'NQR factor present' <==> 'gateway A=7 succeeds'")
else:
    print(f"  COUNTEREXAMPLES FOUND: {a7_with_nqr_fails}")

# === Q4: Does depth >= 2 always guarantee success? ===
print("\n--- Question 4: NQR Depth >= 2 Analysis ---")
print(f"Total Case B QR7 primes analysed: {total_case_b:,}")
print(f"  Saved by depth-1 NQR factor:    {solved_depth1:,}")
print(f"  Saved by depth->=2 NQR factor:  {solved_depth2plus:,}")
if solved_no_nqr_needed:
    print(f"  Anomalies (QA-smooth N saved):  {solved_no_nqr_needed:,}  <-- BUG")

if depth2_saves:
    print(f"\nDepth->=2 saves (first 20):")
    for p, A, d in depth2_saves[:20]:
        N = (p + A) // 4
        factors = factorise(N)
        nqr_info = [(q, e) for q, e in factors.items() if legendre(q, A) == -1]
        print(f"  p={p}, A={A}, N={N}, NQR factors: {nqr_info}")

# Check: for primes saved by depth >= 2, did depth-1 options exist but fail?
print("\n--- Depth-2 saves that had no depth-1 NQR factor anywhere ---")
depth2_no_depth1 = []
for p, A_saver, saved_depth in depth2_saves:
    had_depth1_option = False
    for A in CANDIDATE_AS:
        if (p + A) % 4 != 0:
            continue
        N = (p + A) // 4
        depth, has_nqr = nqr_depth(N, A)
        if has_nqr and depth == 1:
            had_depth1_option = True
            break
    if not had_depth1_option:
        depth2_no_depth1.append((p, A_saver, saved_depth))

print(f"  Cases where depth-2 was essential (no depth-1 anywhere): {len(depth2_no_depth1)}")
for p, A, d in depth2_no_depth1[:10]:
    print(f"    p={p}, A={A}, depth={d}")

# === Per-A statistics ===
print("\n--- Per-A Depth Distribution (top A values) ---")
for A in sorted(a_depth_stats.keys())[:8]:
    stats = a_depth_stats[A]
    total_A = sum(stats.values())
    depth1 = stats.get(1, 0)
    depth2p = sum(v for k, v in stats.items() if k is not None and k >= 2)
    print(f"  A={A:3d}: total={total_A:5d}, depth=1: {depth1:5d} ({100*depth1//max(1,total_A)}%), depth>=2: {depth2p:5d} ({100*depth2p//max(1,total_A)}%)")

print("\nDone.")
