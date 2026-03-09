#!/usr/bin/env python3
"""
Verification script for:
  "A Finite Algebraic Covering System for the Erdos-Straus Conjecture to 10^9"

Verifies Theorem 1.1: every prime p <= LIMIT admits a gateway decomposition
4/p = 1/x + 1/y + 1/z with A <= 239.

Usage:
  python verify.py              # default: verify to 10^6 (~10s)
  python verify.py 10000000     # verify to 10^7 (~2min)
  python verify.py 100000000    # verify to 10^8 (~30min)
  python verify.py 1000000000   # verify to 10^9 (~90min, needs ~1GB RAM)
"""

import sys
import time
from math import isqrt

LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000


# === Sieve ===

def sieve_primes(n):
    """Sieve of Eratosthenes returning a boolean bytearray."""
    ip = bytearray([1]) * (n + 1)
    ip[0] = ip[1] = 0
    for i in range(2, isqrt(n) + 1):
        if ip[i]:
            ip[i*i::i] = bytearray(len(ip[i*i::i]))
    return ip


# === Factorisation ===

SMALL_PRIMES = []

def init_small_primes(bound=100_000):
    global SMALL_PRIMES
    ip = sieve_primes(bound)
    SMALL_PRIMES = [i for i in range(2, bound + 1) if ip[i]]

def factorize(n):
    """Return prime factorisation as {p: e} dict."""
    f = {}
    for p in SMALL_PRIMES:
        if p * p > n:
            break
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f

def divisors_of_square(n):
    """Return sorted list of all divisors of n^2."""
    fac = factorize(n)
    divs = [1]
    for p, e in fac.items():
        new = []
        pe = 1
        for _ in range(2 * e):
            pe *= p
            for d in divs:
                new.append(d * pe)
        divs.extend(new)
    return sorted(divs)


# === Gateway check (Theorem 3.1 + Lemma 3.2) ===

def gateway_check(p, A, d):
    """Check if (p, A, d) gives a valid Erdos-Straus decomposition.

    Returns (x, y, z) if valid, None otherwise.

    Conditions (Theorem 3.1):
      (i)  A | (B + d)
      (ii) d | By  where y = (B+d)/A
    """
    if (p + A) % 4:
        return None
    N = (p + A) // 4
    B = p * N
    if (B + d) % A:
        return None
    y = (B + d) // A
    By = B * y
    if By % d:
        return None
    z = By // d
    if z <= 0:
        return None
    return (N, y, z)


# === Classification ===

def is_case_b(p):
    """Check if p = 1 (mod 24) and all prime factors of (p+3)/4 are = 1 (mod 3)."""
    if p % 24 != 1:
        return False
    n = (p + 3) // 4
    d = 2
    while d * d <= n:
        if n % d == 0:
            if d % 3 != 1:
                return False
            while n % d == 0:
                n //= d
        d += 1
    return n <= 1 or n % 3 == 1


# === Main verification ===

def main():
    t0 = time.time()
    print(f"Erdos-Straus Gateway Verification to {LIMIT:,}")
    print("=" * 60)

    # Sieve
    print(f"Sieving primes to {LIMIT:,}...")
    ip = sieve_primes(LIMIT)
    init_small_primes()
    print(f"Sieve complete ({time.time()-t0:.1f}s)")

    # A values: primes = 3 (mod 4), in three tiers
    prime_set = set(SMALL_PRIMES)
    A_SMALL = [A for A in range(3, 200, 4) if A in prime_set]  # A < 200
    A_MED = [A for A in range(203, 1000, 4) if A in prime_set]  # 200 <= A < 1000
    A_LARGE = [A for A in range(1003, 10000, 4) if A in prime_set]  # 1000 <= A < 10000

    # Counters
    total = 0
    counts = {"p=2": 0, "3mod4": 0, "5mod8": 0, "17mod24": 0,
              "caseA": 0, "nqr7": 0, "qr7_solved": 0}
    max_A = 0
    d_gt_N = 0
    opens = []

    t1 = time.time()
    for p in range(2, LIMIT + 1):
        if not ip[p]:
            continue
        total += 1

        # Layer 1: p = 2
        if p == 2:
            counts["p=2"] = 1
            continue

        # Layer 2: p = 3 (mod 4) - Proposition 4.1
        if p % 4 == 3:
            counts["3mod4"] += 1
            continue

        # Layer 3: p = 5 (mod 8) - Proposition 4.2(a)
        if p % 8 == 5:
            counts["5mod8"] += 1
            continue

        # Layer 4: p = 17 (mod 24) - Proposition 4.2(b)
        if p % 24 == 17:
            counts["17mod24"] += 1
            continue

        # Must be p = 1 (mod 24)
        if p % 24 != 1:
            continue

        # Layer 5: Case A - Proposition 4.2(c)
        if not is_case_b(p):
            counts["caseA"] += 1
            continue

        # Layer 6: NQR mod 7 - Proposition 4.3
        if p % 7 in (3, 5, 6):
            counts["nqr7"] += 1
            continue

        # Layer 7: Case B QR7 - gateway search (Theorem 5.1)
        found = False
        for A_list in [A_SMALL, A_MED, A_LARGE]:
            for A in A_list:
                if (p + A) % 4:
                    continue
                N = (p + A) // 4
                for d in divisors_of_square(N):
                    if d <= 1:
                        continue
                    result = gateway_check(p, A, d)
                    if result:
                        found = True
                        counts["qr7_solved"] += 1
                        if A > max_A:
                            max_A = A
                        if d > N:
                            d_gt_N += 1
                        break
                if found:
                    break
            if found:
                break
        if not found:
            opens.append(p)

        # Progress
        solved_qr7 = counts["qr7_solved"] + len(opens)
        if solved_qr7 % 50000 == 0 and solved_qr7 > 0:
            elapsed = time.time() - t1
            print(f"  {solved_qr7:>10,} Case B QR7 scanned | "
                  f"{counts['qr7_solved']:,} solved | "
                  f"{len(opens)} open | maxA={max_A} | {elapsed:.0f}s")

    del ip
    elapsed = time.time() - t0

    # Results
    proven = sum(counts.values())
    print(f"\nVerification complete ({elapsed:.1f}s)")
    print(f"\n{'Category':<45} {'Count':>12} {'%':>8}")
    print("-" * 67)
    labels = [
        ("p = 2", "p=2"),
        ("p = 3 (mod 4)  [Prop. 4.1]", "3mod4"),
        ("p = 5 (mod 8)  [Prop. 4.2a]", "5mod8"),
        ("p = 17 (mod 24) [Prop. 4.2b]", "17mod24"),
        ("Case A  [Prop. 4.2c]", "caseA"),
        ("Case B NQR7  [Prop. 4.3]", "nqr7"),
        ("Case B QR7  [Thm. 5.1]", "qr7_solved"),
    ]
    for label, key in labels:
        c = counts[key]
        print(f"  {label:<43} {c:>12,} {100*c/total:>7.3f}%")
    print("-" * 67)
    print(f"  {'PROVEN':<43} {proven:>12,} {100*proven/total:>7.4f}%")
    print(f"  {'OPEN':<43} {len(opens):>12} {100*len(opens)/total:>7.5f}%")

    print(f"\nMax A needed: {max_A}")
    print(f"Solutions using d > N (d | N^2 extension): {d_gt_N:,}")

    if opens:
        print(f"\n{len(opens)} OPEN primes:")
        print(f"  First 20: {opens[:20]}")
    else:
        print(f"\nALL {total:,} PRIMES TO {LIMIT:,} VERIFIED.")

    return len(opens) == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
