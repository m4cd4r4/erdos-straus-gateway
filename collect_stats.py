#!/usr/bin/env python3
"""Collect detailed statistics for paper figures at 10^6 scale."""
import json, time
from math import isqrt
from collections import Counter

LIMIT = 1_000_000

def sieve(n):
    ip = bytearray([1]) * (n + 1)
    ip[0] = ip[1] = 0
    for i in range(2, isqrt(n) + 1):
        if ip[i]:
            ip[i*i::i] = bytearray(len(ip[i*i::i]))
    return ip

ip = sieve(LIMIT)
SP = [i for i in range(2, 100000) if ip[i]]
PS = set(SP)

def factorize(n):
    f = {}
    for p in SP:
        if p*p > n: break
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def divs_of_sq(n):
    fac = factorize(n)
    ds = [1]
    for p, e in fac.items():
        new, pe = [], 1
        for _ in range(2*e):
            pe *= p
            for d in ds: new.append(d*pe)
        ds.extend(new)
    return sorted(ds)

def is_case_b(p):
    if p % 24 != 1: return False
    n = (p + 3) // 4
    d = 2
    while d * d <= n:
        if n % d == 0:
            if d % 3 != 1: return False
            while n % d == 0: n //= d
        d += 1
    return n <= 1 or n % 3 == 1

A_ALL = [A for A in range(3, 10000, 4) if A in PS]

# Collect data
stats = {
    "coverage": {},       # layer -> count
    "A_dist": Counter(),  # A -> count (for QR7)
    "d_over_N": [],       # (p, A, d, N, d/N) for each QR7 solution
    "max_A_by_range": {}, # range -> max A seen
}

counts = {"p=2": 0, "3mod4": 0, "5mod8": 0, "17mod24": 0,
          "caseA": 0, "nqr7": 0, "qr7": 0}
total = 0

# Track max A in progressive ranges
range_boundaries = [1000, 2000, 5000, 10000, 20000, 50000,
                    100000, 200000, 500000, 1000000]
max_A_progressive = {}
current_max_A = 0

for p in range(2, LIMIT + 1):
    if not ip[p]: continue
    total += 1
    if p == 2: counts["p=2"] = 1; continue
    if p % 4 == 3: counts["3mod4"] += 1; continue
    if p % 8 == 5: counts["5mod8"] += 1; continue
    if p % 24 == 17: counts["17mod24"] += 1; continue
    if p % 24 != 1: continue
    if not is_case_b(p): counts["caseA"] += 1; continue
    if p % 7 in (3, 5, 6): counts["nqr7"] += 1; continue

    counts["qr7"] += 1
    found = False
    for A in A_ALL:
        if (p + A) % 4: continue
        N = (p + A) // 4
        for d in divs_of_sq(N):
            if d <= 1: continue
            B = p * N
            if (B + d) % A: continue
            y = (B + d) // A
            if (B * y) % d: continue
            z = (B * y) // d
            if z > 0:
                found = True
                stats["A_dist"][A] += 1
                if A > current_max_A:
                    current_max_A = A
                ratio = d / N
                stats["d_over_N"].append({
                    "p": p, "A": A, "d": d, "N": N,
                    "ratio": round(ratio, 4),
                    "d_divides_N": d <= N and N % d == 0
                })
                break
        if found: break
    if not found:
        print(f"OPEN: {p}")

    # Track progressive max A
    for boundary in range_boundaries:
        if p <= boundary and boundary not in max_A_progressive:
            pass  # will set at boundary
    for boundary in range_boundaries:
        if p == boundary or (p > boundary and boundary not in max_A_progressive):
            max_A_progressive[boundary] = current_max_A

# Set any remaining boundaries
for b in range_boundaries:
    if b not in max_A_progressive:
        max_A_progressive[b] = current_max_A

stats["coverage"] = counts
stats["total"] = total
stats["max_A_progressive"] = {str(k): v for k, v in sorted(max_A_progressive.items())}
stats["A_dist"] = {str(k): v for k, v in sorted(stats["A_dist"].items(), key=lambda x: -x[1])}

# Summarise d/N ratios
d_gt_N = sum(1 for x in stats["d_over_N"] if x["ratio"] > 1.0)
stats["d_gt_N_count"] = d_gt_N
stats["d_gt_N_pct"] = round(100 * d_gt_N / max(counts["qr7"], 1), 2)

# Save
with open("stats_1M.json", "w") as f:
    json.dump(stats, f, indent=2, default=str)

print(f"\nCollected stats for {total:,} primes")
print(f"Case B QR7: {counts['qr7']:,}")
print(f"Max A: {current_max_A}")
print(f"d > N cases: {d_gt_N} ({stats['d_gt_N_pct']}%)")
print(f"Saved to stats_1M.json")
