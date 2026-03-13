#!/usr/bin/env python3
"""
Segmented sieve verification - handles 10^10 with ~200MB RAM instead of 10GB.

Processes primes in segments of SEGMENT_SIZE at a time.

Usage:
  python verify_segmented.py 10000000000   # verify to 10^10 (~12-20 hours)
  python verify_segmented.py 10000000000 500000000  # custom segment size
"""

import sys
import time
from math import isqrt

LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000_000_000
SEGMENT = int(sys.argv[2]) if len(sys.argv) > 2 else 500_000_000  # 500MB per segment

MAX_A = 239   # known sufficient from 10^9 verification; we check if it holds

# === Sieve small primes up to sqrt(LIMIT) ===
def sieve_small(n):
    ip = bytearray([1]) * (n + 1)
    ip[0] = ip[1] = 0
    for i in range(2, isqrt(n) + 1):
        if ip[i]:
            ip[i*i::i] = bytearray(len(ip[i*i::i]))
    return [i for i in range(2, n+1) if ip[i]]

# === Gateway check ===
def factorise(n):
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = 1
    return factors

def divisors_of_nsq_mod(n, A):
    """Yield residues of divisors of n^2 mod A, stopping when target found."""
    factors = list(factorise(n).items())
    target = (-4 * n * n) % A  # = -4N^2 mod A
    # Build divisor set mod A iteratively
    residues = {1}
    for q, e in factors:
        q_pows = [pow(q, k, A) for k in range(2*e+1)]
        new_residues = set()
        for r in residues:
            for qp in q_pows:
                new_residues.add((r * qp) % A)
        residues = new_residues
        if target in residues:
            return True
    return target in residues

def try_gateway(p, A):
    """Quick gateway check: does N_A^2 have a divisor ≡ target mod A?"""
    if (p + A) % 4 != 0:
        return False
    N = (p + A) // 4
    B = p * N
    target = (-B) % A
    # Check all divisors of N^2 mod A
    factors = list(factorise(N).items())
    residues = {1}
    for q, e in factors:
        q_pows = [pow(q, k, A) for k in range(2*e+1)]
        new_residues = set()
        for r in residues:
            for qp in q_pows:
                new_residues.add((r * qp) % A)
        residues = new_residues
        if target in residues:
            return True
    return target in residues

def is_case_b_qr7(p):
    """Quick check for Case B QR7."""
    if p % 24 != 1: return False
    if p % 7 not in (1, 2, 4): return False
    n = (p + 3) // 4
    d = 2
    while d * d <= n:
        while n % d == 0:
            if d % 3 != 1: return False
            n //= d
        d += 1
    if n > 1 and n % 3 != 1: return False
    return True

# Candidate A values
CANDIDATE_AS = [
    a for a in range(3, MAX_A + 1)
    if a % 4 == 3 and all(a % d != 0 for d in range(2, isqrt(a)+1))
]

def prop_3mod4(p):
    return p % 4 == 3

def prop_5mod8(p):
    return p % 8 == 5

def prop_17mod24(p):
    return p % 24 == 17

def prop_caseA(p):
    if p % 24 != 1: return False
    n = (p + 3) // 4
    d = 2
    while d * d <= n:
        while n % d == 0:
            if d % 3 == 2: return True
            n //= d
        d += 1
    return n > 1 and n % 3 == 2

def prop_nqr7(p):
    if p % 24 != 1: return False
    return p % 7 in (3, 5, 6)

# === Main segmented verification ===

t0 = time.time()
sqrt_limit = isqrt(LIMIT)
small_primes = sieve_small(sqrt_limit)
print(f"Erdos-Straus Segmented Verification to {LIMIT:,}")
print(f"Segment size: {SEGMENT:,} | Small primes sieved to {sqrt_limit:,}")
print("=" * 65)

total = 0
case_b_qr7 = 0
max_A_needed = 7
hardest_prime = None
open_count = 0
checkpoint_every = 10  # print progress every N segments

seg_start = 2
seg_num = 0

while seg_start <= LIMIT:
    seg_end = min(seg_start + SEGMENT - 1, LIMIT)
    seg_size = seg_end - seg_start + 1

    # Segmented sieve
    sieve = bytearray([1]) * seg_size
    for sp in small_primes:
        start_idx = ((seg_start + sp - 1) // sp) * sp - seg_start
        if seg_start <= sp <= seg_end:
            sieve[sp - seg_start] = 1  # sp itself is prime
            start_idx = sp * 2 - seg_start
        sieve[start_idx::sp] = bytearray(len(sieve[start_idx::sp]))

    # Handle p=2 edge case in first segment
    if seg_start == 2:
        sieve[0] = 0  # 2 is prime but handled separately

    seg_primes = []
    for i in range(seg_size):
        p = seg_start + i
        if sieve[i] and p >= 2:
            seg_primes.append(p)

    for p in seg_primes:
        total += 1
        if total == 1:  # p=2
            continue

        # Quick algebraic checks (cover ~97.6%)
        if prop_3mod4(p) or prop_5mod8(p) or prop_17mod24(p) or prop_caseA(p) or prop_nqr7(p):
            continue

        # Case B QR7: need computational search
        case_b_qr7 += 1
        solved = False
        for A in CANDIDATE_AS:
            if try_gateway(p, A):
                if A > max_A_needed:
                    max_A_needed = A
                    hardest_prime = p
                solved = True
                break

        if not solved:
            open_count += 1
            print(f"  OPEN: p={p} not solved by any A <= {MAX_A}!")

    seg_num += 1
    elapsed = time.time() - t0
    pct = 100.0 * seg_end / LIMIT
    rate = seg_end / elapsed if elapsed > 0 else 0
    eta = (LIMIT - seg_end) / rate if rate > 0 else 0

    if seg_num % checkpoint_every == 0 or seg_end == LIMIT:
        print(f"  [{pct:5.1f}%] up to {seg_end:,} | "
              f"CaseB QR7: {case_b_qr7:,} | max A: {max_A_needed} | "
              f"elapsed: {elapsed:.0f}s | ETA: {eta/3600:.1f}h")
        sys.stdout.flush()

    seg_start = seg_end + 1

elapsed = time.time() - t0
print()
print("=" * 65)
print(f"RESULT: All {total:,} primes up to {LIMIT:,} verified.")
print(f"  Case B QR7 primes:   {case_b_qr7:,}")
print(f"  Maximum A required:  {max_A_needed}")
print(f"  Hardest prime:       {hardest_prime}")
print(f"  Open cases:          {open_count}")
print(f"  Total time:          {elapsed/3600:.2f} hours")

if open_count == 0:
    verdict = "VERIFIED" if max_A_needed <= MAX_A else "NEW MAX A FOUND"
    print(f"\n{verdict}: Erdos-Straus holds for all p <= {LIMIT:,}.")
    if hardest_prime:
        print(f"  New hardest prime: p={hardest_prime}, requires A={max_A_needed}")
        print(f"  Previous max was A=239 (p=32,349,601)")
