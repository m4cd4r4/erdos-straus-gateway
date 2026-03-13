#!/usr/bin/env python3
"""
ERDOS-STRAUS -- SESSION 18: SCALE TO 10,000,000,000 (10^10)
============================================================
Session 17: 100% at 1B (maxA expected ~239-250 based on trend).
Session 18: Push to 10B using GPU-accelerated segmented sieve.

Architecture:
  - Segmented sieve in 10^8 chunks (GPU via CuPy, fallback CPU)
  - Classification: vectorized numpy mod checks + per-prime factorization
  - Gateway scan: CPU with fast trial division for divisors
  - Checkpointing: saves totals after every segment for crash recovery
  - Memory: only one segment in RAM at a time (~100MB)

Key question: Does maxA remain bounded as we push past 1B?
If yes -> finite algebraic covering system through 10^10.
"""
import sys, time, json, os
from math import isqrt
from collections import Counter

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

print("ERDOS-STRAUS -- SESSION 18: SCALE TO 10^10")
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

LIMIT      = 10_000_000_000
SEG        = 100_000_000        # 100M per segment = 100MB sieve
SQRT_LIMIT = isqrt(LIMIT)       # 100,000 exactly
CKPT       = os.path.join(os.path.dirname(__file__), "results", "session18_checkpoint.json")

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
print(f"Small primes to {SQRT_LIMIT:,}: {len(SMALL_PRIMES):,}  ({time.time()-t0:.2f}s)")

# ── A-value lists ────────────────────────────────────────────────────────────

TS0 = [(t, 3+4*t) for t in range(1, 60)   if (3+4*t) in PRIME_SET and (3+4*t) <   200]
TS1 = [(t, 3+4*t) for t in range(1, 260)  if (3+4*t) in PRIME_SET and 200 <= (3+4*t) <  1000]
TS2 = [(t, 3+4*t) for t in range(1, 2600) if (3+4*t) in PRIME_SET and 1000 <= (3+4*t) < 10000]
print(f"A-value phases: {len(TS0)} (<200), {len(TS1)} (<1K), {len(TS2)} (<10K)")

# ── Factorization (trial division, covers N_t up to ~2.5*10^9) ──────────────

def factorize(n):
    if n <= 1: return {}
    f = {}
    for p in SMALL_PRIMES:
        if p * p > n: break
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def divs_of_sq(n):
    fac = factorize(n)
    ds = [1]
    for p, e in fac.items():
        new, pe = [], 1
        for _ in range(2*e):
            pe *= p
            for d in ds:
                new.append(d * pe)
        ds.extend(new)
    return sorted(ds)

def chk(p, t, d):
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

# ── Thm5 case-A check: all prime factors of (p+3)//4 must be ≡ 1 (mod 3) ───

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

# ── Gateway scan for one prime ───────────────────────────────────────────────

def gateway_solve(p):
    """Returns (found, A_used) for a Case B QR7 prime."""
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
    """Returns numpy uint8 array: 1=prime, 0=composite for [lo..hi]."""
    size = hi - lo + 1
    sieve = cp.ones(size, dtype=cp.uint8)
    if lo == 0:
        sieve[0] = 0
    if lo <= 1 <= hi:
        sieve[1 - lo] = 0
    for p in SMALL_PRIMES:
        if p * p > hi:
            break
        start = max(p * p, ((lo + p - 1) // p) * p)
        if start > hi:
            continue
        sieve[start - lo::p] = 0
    return cp.asnumpy(sieve)

def sieve_segment_cpu(lo, hi):
    """CPU fallback segmented sieve."""
    size = hi - lo + 1
    sieve = bytearray([1]) * size
    if lo == 0: sieve[0] = 0
    if lo <= 1 <= hi: sieve[1 - lo] = 0
    for p in SMALL_PRIMES:
        if p * p > hi: break
        start = max(p * p, ((lo + p - 1) // p) * p)
        if start > hi: continue
        for j in range(start - lo, size, p):
            sieve[j] = 0
    return sieve

do_sieve = sieve_segment_gpu if GPU else sieve_segment_cpu

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
    totals = ckpt['totals']
    Au     = Counter({int(k): v for k, v in ckpt['Au'].items()})
    opens  = ckpt['opens']
    mxA    = ckpt['mxA']
else:
    start_lo = 0
    totals   = {'primes': 0, 'p2': 0, 'thm1': 0, 'thm2': 0, 'thm4': 0,
                'thm5': 0, 'thm7': 0, 'gw': 0, 'n_ext': 0}
    Au       = Counter()
    opens    = []
    mxA      = 0

print(f"\nScanning {start_lo:,} → {LIMIT:,} in segments of {SEG:,}")
print("=" * 60)

t_total = time.time()
seg_count = 0

lo = start_lo
while lo <= LIMIT:
    hi = min(lo + SEG - 1, LIMIT)
    t_seg = time.time()

    # ── Sieve this segment ──────────────────────────────────────────────────
    sieve_arr = do_sieve(lo, hi)

    # ── Find primes in segment ──────────────────────────────────────────────
    offsets = np.where(np.frombuffer(sieve_arr, dtype=np.uint8) if isinstance(sieve_arr, bytearray)
                       else sieve_arr)[0]
    primes_in_seg = offsets.astype(np.int64) + lo

    seg_primes = 0
    seg_gw = 0
    seg_ext = 0
    seg_opens = []

    for p in primes_in_seg:
        p = int(p)
        seg_primes += 1
        totals['primes'] += 1

        if p == 2:
            totals['p2'] += 1
            continue

        r4 = p % 4
        if r4 == 3:
            totals['thm1'] += 1
            continue

        if p % 8 == 5:
            totals['thm2'] += 1
            continue

        if p % 24 == 17:
            totals['thm4'] += 1
            continue

        if p % 24 != 1:
            # Should not happen for odd primes — skip (p=2 already handled)
            continue

        # Case A check (Thm 5)
        if is_case_a(p):
            totals['thm5'] += 1
            continue

        # NQR mod 7 check (Thm 7)
        if p % 7 in (3, 5, 6):
            totals['thm7'] += 1
            continue

        # Case B QR7 — needs gateway
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

    seg_time = time.time() - t_seg
    total_solved = totals['gw'] + totals['n_ext']
    total_proven = sum(totals[k] for k in ['p2','thm1','thm2','thm4','thm5','thm7','gw','n_ext'])
    total_open   = len(opens)

    print(f"[{lo//SEG+1:3d}] {lo:>12,}–{hi:>12,}  "
          f"{seg_primes:,} primes  "
          f"gw={seg_gw+seg_ext}  open={len(seg_opens)}  "
          f"maxA={mxA}  {seg_time:.0f}s")
    if seg_opens:
        print(f"  OPEN in this segment: {seg_opens[:5]}")

    # ── Checkpoint after each segment ───────────────────────────────────────
    save_checkpoint({
        'next_lo': hi + 1,
        'totals':  totals,
        'Au':      dict(Au),
        'opens':   opens[:200],   # cap to prevent huge checkpoint files
        'mxA':     mxA,
    })

    lo = hi + 1
    seg_count += 1

# ── Final report ─────────────────────────────────────────────────────────────

elapsed = time.time() - t_total
tot = totals['primes']
total_proven = sum(totals[k] for k in ['p2','thm1','thm2','thm4','thm5','thm7','gw','n_ext'])
total_gw = totals['gw'] + totals['n_ext']

print("\n" + "=" * 60)
print("SCAN COMPLETE")
print("=" * 60)
print(f"  Wall time: {elapsed/3600:.2f}h  ({elapsed:.0f}s)")
print(f"  Segments:  {seg_count}")
print(f"  Total primes to {LIMIT:,}: {tot:,}")

print("\n" + "=" * 60)
print("COVERAGE TABLE (PRIMES TO 10^10)")
print("=" * 60)

labels = [
    ("p=2",                            totals['p2']),
    ("p≡3(mod 4)  [Thm 1]",            totals['thm1']),
    ("p≡5(mod 8)  [Thm 2]",            totals['thm2']),
    ("p≡17(mod 24) [Thm 4]",           totals['thm4']),
    ("Case A      [Thm 5]",            totals['thm5']),
    ("Case B NQR-mod-7 [Thm 7]",       totals['thm7']),
    ("Case B QR-mod-7  [Gateway A<200]", totals['gw']),
    ("Case B QR-mod-7  [Gateway A≥200]", totals['n_ext']),
]
print(f"\n  {'Category':<48}  {'Count':>12}  {'%':>8}")
print("  " + "-"*70)
for lab, v in labels:
    print(f"  {lab:<48}  {v:>12,}  {100*v/max(tot,1):>7.3f}%")
print("  " + "-"*70)
print(f"  {'PROVEN TOTAL':<48}  {total_proven:>12,}  {100*total_proven/max(tot,1):>7.4f}%")
print(f"  {'OPEN':<48}  {len(opens):>12,}  {100*len(opens)/max(tot,1):>7.5f}%")

print("\n" + "=" * 60)
print("A-VALUE DISTRIBUTION (top 30)")
print("=" * 60)
cum = 0
print(f"\n  {'A':>6}  {'count':>10}  {'%':>8}  {'cumul%':>8}")
for A, cnt in Au.most_common(30):
    cum += cnt
    print(f"  {A:>6}  {cnt:>10,}  {100*cnt/max(total_gw,1):>7.2f}%  {100*cum/max(total_gw,1):>7.2f}%")

print("\n" + "=" * 60)
print("SCALE PROGRESSION")
print("=" * 60)
print("  10^6:  78,498 primes       100%  (maxA=199)")
print("  10^7:  664,579 primes      100%  (maxA=199)")
print("  10^8:  5,761,455 primes    100%  (maxA=239)")
print("  10^9:  ~50,847,534 primes  100%  (maxA=?  — see session17 output)")
print(f"  10^10: {tot:,} primes  "
      f"{'100' if not opens else f'{100*total_proven/tot:.4f}'}%  (maxA={mxA})")

if opens:
    print(f"\n{len(opens)} OPEN primes found:")
    print("  First 20:", opens[:20])
else:
    print("\n" + "=" * 60)
    print("ALL PRIMES TO 10,000,000,000 PROVEN ALGEBRAICALLY!")
    print("=" * 60)
