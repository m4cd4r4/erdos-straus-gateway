#!/usr/bin/env python3
"""
ERDŐS-STRAUS -- SESSION 19: SCALE TO 10^11 (exploratory)
=========================================================
Session 18: 100% at 10^10 (maxA=251, n_ext=4, 455,052,511 primes).
Session 19: Push to 10^11 using Pollard-Brent factorization.

Key upgrade over session18:
  Old: O(sqrt(Nt)) trial division through ~316K SMALL_PRIMES
  New: trial division to 1000 (168 primes) + O(Nt^1/4) Pollard-Brent
  Speedup: ~300x per factorize call at the top of the 10^11 range.

NOTE: Exploratory run — results kept separate from the paper.
      Do not update PAPER_PLAN.md or erdos_straus_gateway.tex.
"""
import sys, time, json, os, math
from math import isqrt
from collections import Counter
from random import randint

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

print("ERDŐS-STRAUS -- SESSION 19: SCALE TO 10^11 (exploratory)")
print("=" * 60)

try:
    import cupy as cp
    GPU = True
    print(f"GPU: CuPy {cp.__version__} available")
except ImportError:
    GPU = False
    print("GPU: CuPy not found, using CPU segmented sieve")

import numpy as np

# ── Constants ────────────────────────────────────────────────────────────────

LIMIT      = 100_000_000_000       # 10^11
SEG        = 100_000_000           # 100M per segment
SQRT_LIMIT = isqrt(LIMIT)          # 316,227
CKPT       = os.path.join(os.path.dirname(__file__), "results", "session19_checkpoint.json")

os.makedirs(os.path.dirname(CKPT), exist_ok=True)

# ── Small primes to sqrt(LIMIT) ──────────────────────────────────────────────

t0 = time.time()
ip_small = bytearray([1]) * (SQRT_LIMIT + 1)
ip_small[0] = ip_small[1] = 0
for i in range(2, isqrt(SQRT_LIMIT) + 1):
    if ip_small[i]:
        ip_small[i*i::i] = bytearray(len(ip_small[i*i::i]))
SMALL_PRIMES = [i for i in range(2, SQRT_LIMIT + 1) if ip_small[i]]
PRIME_SET    = set(SMALL_PRIMES)

# Primes used for trial division in factorize (up to 1000 = 168 primes)
TRIAL_PRIMES = [p for p in SMALL_PRIMES if p <= 1000]

print(f"Small primes to {SQRT_LIMIT:,}: {len(SMALL_PRIMES):,}  ({time.time()-t0:.2f}s)")
print(f"Trial-division primes (<=1000): {len(TRIAL_PRIMES)}")

# ── A-value lists ────────────────────────────────────────────────────────────

TS0 = [(t, 3+4*t) for t in range(1, 60)   if (3+4*t) in PRIME_SET and (3+4*t) <   200]
TS1 = [(t, 3+4*t) for t in range(1, 260)  if (3+4*t) in PRIME_SET and 200 <= (3+4*t) <  1000]
TS2 = [(t, 3+4*t) for t in range(1, 2600) if (3+4*t) in PRIME_SET and 1000 <= (3+4*t) < 10000]
print(f"A-value phases: {len(TS0)} (<200), {len(TS1)} (<1K), {len(TS2)} (<10K)")

# ── Fast factorization: trial div to 1000 + Pollard-Brent ───────────────────
#
# For Nt < 2.5e10 (the range we factorize), witnesses {2,3,5,7,11,13} give a
# deterministic Miller-Rabin test (covers all n < 3,474,749,660,383).

def _miller_rabin(n, a):
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return True
    for _ in range(r - 1):
        x = x * x % n
        if x == n - 1:
            return True
    return False

def is_prime_fast(n):
    if n < 2:  return False
    if n < 4:  return True
    if n % 2 == 0: return False
    for a in (2, 3, 5, 7, 11, 13):
        if n == a: return True
        if not _miller_rabin(n, a): return False
    return True

def _brent_rho(n):
    """Brent's rho: returns a non-trivial factor of composite n.

    Backtracking loop is bounded to prevent infinite spin when ys hits a
    fixed point of f(x)=x²+c mod n (which gives gcd(0,n)=n forever).
    """
    if n % 2 == 0:
        return 2
    while True:
        y = randint(1, n - 1)
        c = randint(1, n - 1)
        m = randint(1, n - 1)
        g = q = r = 1
        x = ys = 0
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r - k)):
                    y  = (y * y + c) % n
                    q  = q * abs(x - y) % n
                g  = math.gcd(q, n)
                k += m
            r *= 2
        if g == n:
            # Backtrack step by step; bound the loop to avoid fixed-point spin
            g = n
            steps = 0
            while g == n and steps < 1000:
                ys = (ys * ys + c) % n
                g  = math.gcd(abs(x - ys), n)
                steps += 1
        if g != n:
            return g
        # unlucky c or hit fixed point — retry with new random values

def _factor_large(n, f):
    """Recursively fully factor n into primes, accumulating in dict f."""
    if n <= 1:
        return
    if is_prime_fast(n):
        f[n] = f.get(n, 0) + 1
        return
    d = _brent_rho(n)
    _factor_large(d, f)
    _factor_large(n // d, f)

def factorize(n):
    if n <= 1:
        return {}
    f = {}
    for p in TRIAL_PRIMES:
        if p * p > n:
            break
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
    if n > 1:
        _factor_large(n, f)
    return f

def divs_of_sq(n):
    fac = factorize(n)
    ds = [1]
    for p, e in fac.items():
        new, pe = [], 1
        for _ in range(2 * e):
            pe *= p
            for d in ds:
                new.append(d * pe)
        ds.extend(new)
    return sorted(ds)

# ── Gateway check ────────────────────────────────────────────────────────────

def chk(p, t, d):
    A  = 3 + 4 * t
    xn = p + A
    if xn % 4: return False
    x  = xn // 4
    B  = p * x
    if (B + d) % A: return False
    y  = (B + d) // A
    By = B * y
    if By % d: return False
    z  = By // d
    return z > 0 and 4 * x * y * z == p * (y * z + x * z + x * y)

def is_case_a(p):
    n2 = (p + 3) // 4
    d2 = 2
    while d2 * d2 <= n2:
        if n2 % d2 == 0:
            if d2 % 3 != 1:
                return False
            while n2 % d2 == 0:
                n2 //= d2
        d2 += 1
    return not (n2 > 1 and n2 % 3 != 1)

def gateway_solve(p):
    for phase in [TS0, TS1, TS2]:
        for (t, A) in phase:
            if (p + A) % 4: continue
            Nt = (p + A) // 4
            for d in divs_of_sq(Nt):
                if d <= 1: continue
                if chk(p, t, d):
                    return True, A
    return False, 0

# ── GPU segmented sieve ──────────────────────────────────────────────────────

def sieve_segment_gpu(lo, hi):
    size  = hi - lo + 1
    sieve = cp.ones(size, dtype=cp.uint8)
    if lo == 0:        sieve[0] = 0
    if lo <= 1 <= hi:  sieve[1 - lo] = 0
    for p in SMALL_PRIMES:
        if p * p > hi: break
        start = max(p * p, ((lo + p - 1) // p) * p)
        if start > hi: continue
        sieve[start - lo::p] = 0
    return cp.asnumpy(sieve)

def sieve_segment_cpu(lo, hi):
    size  = hi - lo + 1
    sieve = bytearray([1]) * size
    if lo == 0:        sieve[0] = 0
    if lo <= 1 <= hi:  sieve[1 - lo] = 0
    for p in SMALL_PRIMES:
        if p * p > hi: break
        start = max(p * p, ((lo + p - 1) // p) * p)
        if start > hi: continue
        for j in range(start - lo, size, p):
            sieve[j] = 0
    return sieve

do_sieve = sieve_segment_gpu if GPU else sieve_segment_cpu

# ── Quick validation (10M) ───────────────────────────────────────────────────

print("\nValidating factorization at 10M...")
val_sieve = sieve_segment_cpu(0, 10_000_000)
val_opens = []
for off in range(len(val_sieve)):
    if not val_sieve[off]: continue
    p = int(off)
    if p < 2: continue
    if p == 2 or p % 4 == 3 or p % 8 == 5 or p % 24 == 17: continue
    if p % 24 != 1: continue
    if is_case_a(p): continue
    if p % 7 in (3, 5, 6): continue
    found, _ = gateway_solve(p)
    if not found:
        val_opens.append(p)

if val_opens:
    print(f"  VALIDATION FAILED: {len(val_opens)} opens at 10M: {val_opens[:5]}")
    sys.exit(1)
else:
    print("  OK: 0 opens at 10M. Fast factorization matches session18.")

# ── Checkpoint helpers ───────────────────────────────────────────────────────

def load_checkpoint():
    if os.path.exists(CKPT):
        with open(CKPT) as f:
            return json.load(f)
    return None

def save_checkpoint(state):
    with open(CKPT, 'w') as f:
        json.dump(state, f)

# ── Main scan ────────────────────────────────────────────────────────────────

ckpt = load_checkpoint()
if ckpt:
    print(f"\nRESUMING from checkpoint: segment starting at {ckpt['next_lo']:,}")
    start_lo = ckpt['next_lo']
    totals   = ckpt['totals']
    Au       = Counter({int(k): v for k, v in ckpt['Au'].items()})
    opens    = ckpt['opens']
    mxA      = ckpt['mxA']
else:
    start_lo = 0
    totals   = {'primes': 0, 'p2': 0, 'thm1': 0, 'thm2': 0, 'thm4': 0,
                'thm5': 0, 'thm7': 0, 'gw': 0, 'n_ext': 0}
    Au       = Counter()
    opens    = []
    mxA      = 0

print(f"\nScanning {start_lo:,} -> {LIMIT:,} in segments of {SEG:,}  ({(LIMIT-start_lo)//SEG} segments)")
print("=" * 60)

t_total = time.time()

lo = start_lo
while lo <= LIMIT:
    hi      = min(lo + SEG - 1, LIMIT)
    t_seg   = time.time()

    sieve_arr = do_sieve(lo, hi)

    offsets       = np.where(np.frombuffer(sieve_arr, dtype=np.uint8)
                             if isinstance(sieve_arr, bytearray) else sieve_arr)[0]
    primes_in_seg = offsets.astype(np.int64) + lo

    seg_primes = seg_gw = seg_ext = 0
    seg_opens  = []

    for p in primes_in_seg:
        p = int(p)
        seg_primes      += 1
        totals['primes'] += 1

        if p == 2:
            totals['p2'] += 1
            continue

        if p % 4 == 3:
            totals['thm1'] += 1
            continue

        if p % 8 == 5:
            totals['thm2'] += 1
            continue

        if p % 24 == 17:
            totals['thm4'] += 1
            continue

        if p % 24 != 1:
            continue

        if is_case_a(p):
            totals['thm5'] += 1
            continue

        if p % 7 in (3, 5, 6):
            totals['thm7'] += 1
            continue

        found, A_used = gateway_solve(p)
        if found:
            Au[A_used] += 1
            if A_used > mxA:
                mxA = A_used
            if A_used < 200:
                totals['gw'] += 1
                seg_gw += 1
            else:
                totals['n_ext'] += 1
                seg_ext += 1
        else:
            opens.append(p)
            seg_opens.append(p)

    seg_time    = time.time() - t_seg
    total_proven = sum(totals[k] for k in ['p2','thm1','thm2','thm4','thm5','thm7','gw','n_ext'])

    print(f"[{lo//SEG+1:4d}] {lo:>14,}-{hi:>14,}  "
          f"{seg_primes:,} primes  "
          f"gw={seg_gw+seg_ext}  open={len(seg_opens)}  "
          f"maxA={mxA}  {seg_time:.0f}s")
    if seg_opens:
        print(f"  OPEN in this segment: {seg_opens[:5]}")

    save_checkpoint({
        'next_lo': hi + 1,
        'totals':  totals,
        'Au':      dict(Au),
        'opens':   opens[:200],
        'mxA':     mxA,
    })

    lo = hi + 1

# ── Final report ─────────────────────────────────────────────────────────────

elapsed = time.time() - t_total
tot     = totals['primes']
total_proven = sum(totals[k] for k in ['p2','thm1','thm2','thm4','thm5','thm7','gw','n_ext'])
total_gw     = totals['gw'] + totals['n_ext']

print("\n" + "=" * 60)
print("SCAN COMPLETE")
print("=" * 60)
print(f"  Wall time : {elapsed/3600:.2f}h  ({elapsed:.0f}s)")
print(f"  Segments  : {(LIMIT - start_lo) // SEG}")
print(f"  Total primes to {LIMIT:,}: {tot:,}")

print("\n" + "=" * 60)
print("COVERAGE TABLE (PRIMES TO 10^11)")
print("=" * 60)
labels = [
    ("p=2",                              totals['p2']),
    ("p≡3(mod 4)  [Thm 1]",             totals['thm1']),
    ("p≡5(mod 8)  [Thm 2]",             totals['thm2']),
    ("p≡17(mod 24) [Thm 4]",            totals['thm4']),
    ("Case A      [Thm 5]",             totals['thm5']),
    ("Case B NQR-mod-7 [Thm 7]",        totals['thm7']),
    ("Case B QR-mod-7  [Gateway A<200]", totals['gw']),
    ("Case B QR-mod-7  [Gateway A>=200]",totals['n_ext']),
]
print(f"\n  {'Category':<48}  {'Count':>12}  {'%':>8}")
print("  " + "-" * 70)
for lab, v in labels:
    print(f"  {lab:<48}  {v:>12,}  {100*v/max(tot,1):>7.3f}%")
print("  " + "-" * 70)
print(f"  {'PROVEN TOTAL':<48}  {total_proven:>12,}  {100*total_proven/max(tot,1):>7.4f}%")
print(f"  {'OPEN':<48}  {len(opens):>12,}  {100*len(opens)/max(tot,1):>7.5f}%")

print("\n" + "=" * 60)
print("A-VALUE DISTRIBUTION (top 30)")
print("=" * 60)
cum = 0
print(f"\n  {'A':>6}  {'count':>12}  {'%':>8}  {'cumul%':>8}")
for A, cnt in Au.most_common(30):
    cum += cnt
    print(f"  {A:>6}  {cnt:>12,}  {100*cnt/max(total_gw,1):>7.2f}%  {100*cum/max(total_gw,1):>7.2f}%")

print("\n" + "=" * 60)
print("SCALE PROGRESSION")
print("=" * 60)
print("  10^6 :       78,498 primes  100%  maxA=199")
print("  10^7 :      664,579 primes  100%  maxA=199")
print("  10^8 :    5,761,455 primes  100%  maxA=239")
print("  10^9 :   50,847,534 primes  100%  maxA=239")
print("  10^10:  455,052,511 primes  100%  maxA=251  n_ext=4")
pct = "100" if not opens else f"{100*total_proven/tot:.4f}"
print(f"  10^11: {tot:,} primes  {pct}%  maxA={mxA}  n_ext={totals['n_ext']}")

if opens:
    print(f"\n{len(opens)} OPEN primes:")
    print("  First 20:", opens[:20])
else:
    print("\n" + "=" * 60)
    print("ALL PRIMES TO 100,000,000,000 PROVEN ALGEBRAICALLY!")
    print("=" * 60)
