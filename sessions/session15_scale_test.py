#!/usr/bin/env python3
"""
ERDOS-STRAUS -- SESSION 15: SCALE TEST TO 10,000,000
=====================================================

Session 14 achieved 100% algebraic coverage for primes to 1,000,000
using Theorem 14 (d | N_t^2 gateway). Session 15 tests whether
this framework scales to 10M.

Questions:
  1. Does the corrected gateway (d | N_t^2, A < 200) cover all
     Case B QR-mod-7 primes to 10M?
  2. If not, how many remain open, and what A_t range closes them?
  3. Theoretical: can we bound the required A_t as a function of p?

Structure:
  Section 1: Full classification of primes to 10M
  Section 2: Gateway scan (d | N_t^2, A < 200) for Case B QR7 primes
  Section 3: Extended A_t for any remaining open primes
  Section 4: Coverage table + analysis
  Section 5: Theoretical observations toward a full proof
"""

import sys
import time
from math import isqrt
from collections import Counter

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

print("ERDOS-STRAUS -- SESSION 15: SCALE TEST TO 10,000,000")
print("=" * 60)

# ── Sieves ───────────────────────────────────────────────────────────────────

LIMIT = 10_000_000

def sieve(n):
    ip = bytearray([1]) * (n+1); ip[0] = ip[1] = 0
    for i in range(2, isqrt(n)+1):
        if ip[i]: ip[i*i::i] = bytearray(len(ip[i*i::i]))
    return [i for i in range(2, n+1) if ip[i]]

t0 = time.time()
PRIMES = sieve(LIMIT)
PRIME_SET = set(PRIMES)
print(f"Sieve to {LIMIT:,}: {len(PRIMES):,} primes ({time.time()-t0:.1f}s)")

# SPF sieve for factorizing N_t values (up to ~2.6M)
SPF_LIMIT = 3_000_001
spf = list(range(SPF_LIMIT))
for i in range(2, isqrt(SPF_LIMIT)+1):
    if spf[i] == i:
        for j in range(i*i, SPF_LIMIT, i):
            if spf[j] == j: spf[j] = i

def factorize(n):
    if n <= 1: return {}
    if n < SPF_LIMIT:
        f = {}
        while n > 1:
            p = spf[n]; f[p] = f.get(p,0)+1; n //= p
        return f
    f = {}; d = 2
    while d*d <= n:
        while n % d == 0: f[d] = f.get(d,0)+1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n,0)+1
    return f

def all_divisors_of_square(n):
    fac = factorize(n)
    divs = [1]
    for p, e in fac.items():
        new = []
        pe = 1
        for k in range(1, 2*e+1):
            pe *= p
            for d in divs:
                new.append(d * pe)
        divs.extend(new)
    return sorted(divs)

# ── Core formula ──────────────────────────────────────────────────────────────

def check_formula_v2(p, t, d):
    """Corrected gateway: d | N_t^2 (not d | N_t)."""
    A  = 3 + 4*t
    xn = p + A
    if xn % 4: return False
    x  = xn // 4
    B  = p * x
    if (B + d) % A: return False
    y  = (B + d) // A
    By = B * y
    if By % d: return False
    z  = By // d
    return z > 0 and 4*x*y*z == p*(y*z + x*z + x*y)

# Prime A_t values
TS_200  = [t for t in range(1, 60) if (3+4*t) in PRIME_SET and (3+4*t) < 200]
TS_1000 = [t for t in range(1, 260) if (3+4*t) in PRIME_SET and (3+4*t) < 1000]
TS_10K  = [t for t in range(1, 2600) if (3+4*t) in PRIME_SET and (3+4*t) < 10000]

# ── Section 1: Classify all primes ───────────────────────────────────────────

print("\n" + "=" * 60)
print("SECTION 1: CLASSIFY PRIMES TO 10M")
print("=" * 60)

t0 = time.time()

# Categories
n_p2    = 1   # p=2
n_thm1  = 0   # p=3 mod 4
n_thm2  = 0   # p=5 mod 8
n_thm4  = 0   # p=17 mod 24
n_thm5  = 0   # Case A
n_thm7  = 0   # Case B NQR7
case_b_qr7 = []

for p in PRIMES:
    if p == 2: continue
    if p % 4 == 3:
        n_thm1 += 1; continue
    if p % 8 == 5:
        n_thm2 += 1; continue
    if p % 24 == 17:
        n_thm4 += 1; continue
    # p = 1 mod 24
    if p % 24 != 1: continue  # shouldn't happen
    tmp = (p + 3) // 4
    all_1mod3 = True
    d2 = 2
    n2 = tmp
    while d2*d2 <= n2:
        if n2 % d2 == 0:
            if d2 % 3 != 1:
                all_1mod3 = False; break
            while n2 % d2 == 0: n2 //= d2
        d2 += 1
    if all_1mod3 and n2 > 1 and n2 % 3 != 1:
        all_1mod3 = False

    if not all_1mod3:
        n_thm5 += 1  # Case A
    elif p % 7 in (3, 5, 6):
        n_thm7 += 1  # Case B NQR7
    else:
        case_b_qr7.append(p)  # Case B QR7 — the hard primes

print(f"Classification complete ({time.time()-t0:.1f}s)")
print(f"  p=2:           {n_p2}")
print(f"  Thm 1 (3mod4): {n_thm1}")
print(f"  Thm 2 (5mod8): {n_thm2}")
print(f"  Thm 4 (17m24): {n_thm4}")
print(f"  Case A (Thm5): {n_thm5}")
print(f"  CaseB NQR7 (7):{n_thm7}")
print(f"  CaseB QR7:     {len(case_b_qr7)}  <-- target for gateway scan")
total = len(PRIMES)
accounted = n_p2 + n_thm1 + n_thm2 + n_thm4 + n_thm5 + n_thm7 + len(case_b_qr7)
print(f"  Total:         {total} (accounted: {accounted}, match: {total == accounted})")

# ── Section 2: Gateway scan (A < 200) ────────────────────────────────────────

print("\n" + "=" * 60)
print("SECTION 2: GATEWAY SCAN (d | N_t^2, A_t < 200)")
print("=" * 60)

t0 = time.time()
open_after_200 = []
solved_200 = 0

for idx, p in enumerate(case_b_qr7):
    covered = False
    for t in TS_200:
        A = 3 + 4*t
        if (p + A) % 4: continue
        Nt = (p + A) // 4
        for d in all_divisors_of_square(Nt):
            if d <= 1: continue
            if check_formula_v2(p, t, d):
                covered = True; break
        if covered: break
    if covered:
        solved_200 += 1
    else:
        open_after_200.append(p)

    if (idx+1) % 5000 == 0:
        elapsed = time.time() - t0
        print(f"  ... {idx+1}/{len(case_b_qr7)} scanned, "
              f"{solved_200} solved, {len(open_after_200)} open ({elapsed:.1f}s)")

elapsed = time.time() - t0
print(f"\nGateway scan (A < 200) complete ({elapsed:.1f}s)")
print(f"  Solved: {solved_200} / {len(case_b_qr7)}")
print(f"  Open:   {len(open_after_200)}")
if open_after_200:
    print(f"  First 30 open: {open_after_200[:30]}")

# ── Section 3: Extended A_t for open primes ──────────────────────────────────

if open_after_200:
    print("\n" + "=" * 60)
    print("SECTION 3: EXTENDED A_t SEARCH FOR OPEN PRIMES")
    print("=" * 60)

    # Phase 1: A < 1000
    open_after_1K = []
    solved_1K = 0
    t0 = time.time()
    for p in open_after_200:
        covered = False
        for t in TS_1000:
            A = 3+4*t
            if A < 200: continue
            if (p+A) % 4: continue
            Nt = (p+A)//4
            for d in all_divisors_of_square(Nt):
                if d <= 1: continue
                if check_formula_v2(p, t, d):
                    covered = True; break
            if covered: break
        if covered:
            solved_1K += 1
        else:
            open_after_1K.append(p)

    print(f"Phase 1 (A < 1000): solved {solved_1K}, open {len(open_after_1K)} ({time.time()-t0:.1f}s)")

    # Phase 2: A < 10000
    if open_after_1K:
        open_after_10K = []
        solved_10K = 0
        t0 = time.time()
        for p in open_after_1K:
            covered = False
            for t in TS_10K:
                A = 3+4*t
                if A < 1000: continue
                if (p+A) % 4: continue
                Nt = (p+A)//4
                for d in all_divisors_of_square(Nt):
                    if d <= 1: continue
                    if check_formula_v2(p, t, d):
                        covered = True; break
                if covered: break
            if covered:
                solved_10K += 1
            else:
                open_after_10K.append(p)

        print(f"Phase 2 (A < 10000): solved {solved_10K}, open {len(open_after_10K)} ({time.time()-t0:.1f}s)")
    else:
        open_after_10K = []

    # Phase 3: A < 100000
    if open_after_10K:
        TS_100K = [t for t in range(1, 25001) if (3+4*t) in PRIME_SET and (3+4*t) < 100000]
        open_after_100K = []
        solved_100K = 0
        t0 = time.time()
        for p in open_after_10K:
            covered = False
            for t in TS_100K:
                A = 3+4*t
                if A < 10000: continue
                if (p+A) % 4: continue
                Nt = (p+A)//4
                for d in all_divisors_of_square(Nt):
                    if d <= 1: continue
                    if check_formula_v2(p, t, d):
                        covered = True; break
                if covered: break
            if covered:
                solved_100K += 1
            else:
                open_after_100K.append(p)

        print(f"Phase 3 (A < 100000): solved {solved_100K}, open {len(open_after_100K)} ({time.time()-t0:.1f}s)")
    else:
        open_after_100K = []

    final_open = open_after_100K if open_after_10K else (open_after_1K if not open_after_1K else [])
    if not open_after_200:
        final_open = []
    elif not open_after_1K:
        final_open = []
    elif not open_after_10K:
        final_open = []
    else:
        final_open = open_after_100K
else:
    final_open = []

# ── Section 4: Coverage table ─────────────────────────────────────────────────

print("\n" + "=" * 60)
print("SECTION 4: FINAL COVERAGE TABLE (PRIMES TO 10M)")
print("=" * 60)

n_gateway = len(case_b_qr7) - len(final_open)
proven = n_p2 + n_thm1 + n_thm2 + n_thm4 + n_thm5 + n_thm7 + n_gateway

print(f"\n  {'Category':<48}  {'Count':>9}  {'%':>8}  Status")
print(f"  {'-'*48}  {'-'*9}  {'-'*8}  ----------")
print(f"  {'p=2':<48}  {n_p2:>9}  {100*n_p2/total:>7.3f}%  PROVEN")
print(f"  {'p=3(mod 4) [Thm 1]':<48}  {n_thm1:>9}  {100*n_thm1/total:>7.3f}%  PROVEN")
print(f"  {'p=5(mod 8) [Thm 2]':<48}  {n_thm2:>9}  {100*n_thm2/total:>7.3f}%  PROVEN")
print(f"  {'p=17(mod 24) [Thm 4]':<48}  {n_thm4:>9}  {100*n_thm4/total:>7.3f}%  PROVEN")
print(f"  {'Case A [Thm 5]':<48}  {n_thm5:>9}  {100*n_thm5/total:>7.3f}%  PROVEN")
print(f"  {'Case B NQR-mod-7 [Thm 7]':<48}  {n_thm7:>9}  {100*n_thm7/total:>7.3f}%  PROVEN")
print(f"  {'Case B QR-mod-7 [Thms 9-14]':<48}  {n_gateway:>9}  {100*n_gateway/total:>7.3f}%  PROVEN")
print(f"  {'-'*48}  {'-'*9}  {'-'*8}  ----------")
print(f"  {'PROVEN TOTAL':<48}  {proven:>9}  {100*proven/total:>7.4f}%")
print(f"  {'OPEN':<48}  {len(final_open):>9}  {100*len(final_open)/total:>7.4f}%")
print(f"  {'-'*48}  {'-'*9}  {'-'*8}  ----------")
print(f"  {'ALL PRIMES':<48}  {total:>9}  {'100.000%':>8}")

# ── Section 5: Analysis of open primes (if any) ─────────────────────────────

if final_open:
    print("\n" + "=" * 60)
    print("SECTION 5: OPEN PRIMES ANALYSIS")
    print("=" * 60)
    print(f"\n{len(final_open)} primes remain open after A_t < 100,000.")
    print(f"First 20: {final_open[:20]}")

    # Factorize N_t for t=1 (A=7)
    print("\nN_t=(p+7)/4 factorizations for first 15 open primes:")
    for p in final_open[:15]:
        Nt = (p+7)//4
        fac = factorize(Nt)
        fstr = " * ".join(f"{q}^{e}" if e>1 else str(q) for q,e in sorted(fac.items()))
        print(f"  p={p:>10}, N_t={Nt:>10} = {fstr}")

# ── Section 6: Theoretical observations ──────────────────────────────────────

print("\n" + "=" * 60)
print("SECTION 6: THEORETICAL OBSERVATIONS")
print("=" * 60)

print(f"""
COMPUTATIONAL EVIDENCE:
  Primes to  1,000,000: 100.000% proven (Session 14)
  Primes to 10,000,000: {100*proven/total:.4f}% proven (Session 15)
  Open: {len(final_open)} primes

THEOREM 14 (Corrected Gateway):
  For prime A=3+4t, d with gcd(d,pA)=1, the formula works when:
    (i)  4 | p+A
    (ii) A | p*N_t + d    [where N_t = (p+A)/4]
    (iii) d | N_t^2

  Equivalently: d | N_t^2 AND d = -p^2/4 (mod A).

PATH TO FULL PROOF:
  Need to show: for every prime p=1(mod 24) with p QR mod 7 [Case B QR7],
  there exists a prime A=3+4t and d | N_t^2 satisfying condition (ii).

  Key observations:
  1. N_t = (p+A)/4 varies with A, giving different factorizations
  2. The divisors of N_t^2 are richer than N_t (higher prime powers)
  3. As A grows, N_t grows, gaining more prime factors
  4. By Chebotarev/Dirichlet, for each A, the residue class d=-p^2/4 (mod A)
     is "hit" by divisors of N_t^2 with positive probability
  5. The probability of ESCAPING all A_t < B decays as B grows

  CONJECTURE (Strong Form):
    For every prime p > 2, there exists a prime A < C*log(p)^2
    and d | ((p+A)/4)^2 such that the gateway formula gives
    4/p = 1/x + 1/y + 1/z with x,y,z positive integers.

  This would imply the Erdos-Straus conjecture.
""")

if not final_open:
    print("=" * 60)
    print("ALL PRIMES TO 10,000,000 PROVEN ALGEBRAICALLY!")
    print("=" * 60)
