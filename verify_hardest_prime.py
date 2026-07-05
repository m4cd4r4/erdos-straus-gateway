#!/usr/bin/env python3
"""
Reproduce every claim in the paper's "Anatomy of the hardest prime" section
for p = 3,807,728,761 (the unique prime < 10^11 requiring A = 359).

Standard library only. Run:  python verify_hardest_prime.py
"""

def factorize(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f

def divisor_residues(n, m):
    """{ d mod m : d divides n^2 }, excluding 0."""
    res = {1}
    for q, e in factorize(n).items():
        nxt, qp = set(), 1
        for _ in range(2 * e + 1):
            for r in res:
                nxt.add((r * qp) % m)
            qp = (qp * q) % m
        res = nxt
    return res - {0}

def order(a, m):
    a %= m
    x, k = a, 1
    while x != 1:
        x = (x * a) % m
        k += 1
    return k

def check(label, got, want):
    ok = got == want
    print(f"  [{'OK' if ok else 'FAIL'}] {label}: {got}" + ("" if ok else f"  (expected {want})"))
    assert ok, label

p, A = 3_807_728_761, 359
N = (p + A) // 4
print(f"p = {p:,}   A = {A}   N = (p+A)/4 = {N:,}")

print("\nWhy p is Case B:")
check("(p+3)/4 factorisation", factorize((p + 3) // 4), {7: 1, 73: 1, 1_862_881: 1})
check("all factors = 1 (mod 3)", all(q % 3 == 1 for q in factorize((p + 3) // 4)), True)

print("\nWhy every A < 359 fails:")
def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n**0.5) + 1))
cands = [a for a in range(7, 359, 4) if is_prime(a)]
worst = 0
for a in cands:
    Na = (p + a) // 4
    t = (-pow(p, 2, a) * pow(4, -1, a)) % a
    R = divisor_residues(Na, a)
    assert t not in R, f"A={a} unexpectedly succeeds"
    worst = max(worst, len(R))
check("candidate A in [7,359) all fail", len(cands), 35)
check("max divisor-residue coverage |R_A|", worst, 184)

print("\nWhy A = 359 succeeds:")
check("N_359 factorisation", factorize(N), {2: 3, 3: 1, 5: 1, 13: 1, 23: 1, 43: 1, 617: 1})
S = divisor_residues(N, A)
check("|R_359| = all 358 nonzero residues", len(S), 358)
t = (-pow(p, 2, A) * pow(4, -1, A)) % A
check("target residue t", t, 140)
d = 1935
check("d = 1935 hits the target", d % A, t)
B = p * N; y = (B + d) // A; z = (B * y) // d
check("identity 4/p = 1/x+1/y+1/z", 4 * N * y * z, p * (y * z + N * z + N * y))
print(f"       x = {N:,}\n       y = {y:,}\n       z = {z:,}")

print("\nThe d | N^2, d not | N phenomenon:")
check("d divides N^2", (N * N) % d, 0)
check("d does not divide N", N % d != 0, True)

print("\nWhy full coverage occurs (discrete-log structure, base g=13):")
check("13 is a primitive root mod 359", order(13, A), 358)
check("{2,3,5,23} are quadratic residues (even index)",
      all(order(q, A) == 179 for q in (2, 3, 5, 23)), True)
check("{13,43,617} are primitive roots (odd index)",
      all(order(q % A, A) == 358 for q in (13, 43, 617)), True)

print("\nSecond mechanism (A=311, p=6,372,847,849):")
p2, A2 = 6_372_847_849, 311
N2 = (p2 + A2) // 4
check("N_311 factorisation", factorize(N2), {2: 3, 3: 2, 5: 1, 7: 1, 632_227: 1})
check("632227 has order 10 (not a primitive root)", order(632_227, A2), 10)
check("coverage without 632227", len(divisor_residues(N2 // 632_227, A2)), 155)
check("coverage with 632227 (full)", len(divisor_residues(N2, A2)), 310)

print("\nAll anatomy claims verified.")
