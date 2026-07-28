#!/usr/bin/env python3
"""
ERDOS-STRAUS -- SESSION 22: POLLARD-BRENT FACTORISATION ENGINE
================================================================
Ports session19_100B.py's fast factoriser (trial division to 1000, then
Miller-Rabin + Pollard-Brent for the residual cofactor) onto session20_
corrected_10B.py's engine, WITHOUT touching its classify()/gateway logic.

Why: at 10^11, N=(p+A)/4 is up to ~2.5e10 and trial division to sqrt(N)
(~160K primes) is affordable. At 10^12, N is up to ~2.5e11 and sqrt(N)
~500K primes per hard factorisation -- the classifier itself (the Case
A/B check on (p+3)/4) and the gateway divisor search (on N^2) both pay
this cost per prime. Trial division to 1000 (168 primes) peels off small
factors cheaply; Miller-Rabin certifies primality of whatever remains;
Pollard-Brent (O(N^1/4)) finds a factor of any residual composite. This
also speeds up the Case A/B check itself, which session19 never did
(it kept plain trial division there) -- both hot paths get the same
treatment here.

Correctness invariant: classify()'s MATH is byte-identical to session20's.
Only the factorisation backend changed. Since _divs_of_sq always returns
a sorted() divisor list and the gateway search takes the first divisor
(in ascending numeric order) satisfying chk(), the (p, A, d) triple
returned for any solved prime is uniquely determined by the true
divisor set of N^2 -- independent of which algorithm found the factors.
So this engine must reproduce session20/21's checkpoints EXACTLY
(same counts, same A-distribution, same witnesses) at every milestone
below 10^11. See results/session22_validation_checkpoint.json runs.

Miller-Rabin witness bound: {2,3,5,7,11,13} is a PROVEN deterministic
test for n < 3,474,749,660,383 (~3.47e12) [Jaeschke 1993]. At LIMIT=10^12,
N=(p+A)/4 <= ~2.5e11, safely under this bound. Do not push this engine
past ~3e12 without extending the witness set (add 17 for deterministic
coverage to 3.41e14).

Usage:
  python session22_pollard_brent.py [LIMIT] [CKPT_PATH]
  Defaults: LIMIT = 10^10 (safe default; pass 10^12 explicitly for the
  real attempt), CKPT = results/session22_checkpoint.json.
"""
import sys, time, json, os, math
from math import isqrt
from collections import Counter
from multiprocessing import Pool
from random import randint

LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000_000_000
CKPT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "results", "session22_checkpoint.json")

SEG        = 100_000_000
SQRT_LIMIT = max(isqrt(LIMIT), 1000)   # sieve bound only; factoriser no longer needs it
WORKERS    = int(os.environ.get("ES_WORKERS", 16))
# Half-decade milestones matching session20_corrected_10B.py (2026-07-25/26).
MILESTONES = [10**6, 10**7, 10**8, 10**9, 10**10, 10**11,
              2 * 10**11, 5 * 10**11, 10**12,
              2 * 10**12, 5 * 10**12, 10**13]

# ── Worker globals (built once per spawned process via initializer) ─────────

_TRIAL = None   # primes <= 1000, for cheap trial division (168 primes)
_TS = None      # candidate A values: primes = 3 (mod 4), ascending, < 10^4

def _small_primes(n):
    ip = bytearray([1]) * (n + 1)
    ip[0] = ip[1] = 0
    for i in range(2, isqrt(n) + 1):
        if ip[i]:
            ip[i*i::i] = bytearray(len(ip[i*i::i]))
    return [i for i in range(2, n + 1) if ip[i]]

def _init_worker():
    global _TRIAL, _TS
    _TRIAL = _small_primes(1000)
    _TS = [A for A in _small_primes(10_000) if A % 4 == 3]

# ── Fast factorisation: trial div to 1000 + Miller-Rabin + Pollard-Brent ────

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

def _is_prime_fast(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    for a in (2, 3, 5, 7, 11, 13):
        if n == a: return True
        if not _miller_rabin(n, a): return False
    return True

def _brent_rho(n):
    """Returns a non-trivial factor of composite n (Brent's rho variant)."""
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
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = math.gcd(q, n)
                k += m
            r *= 2
        if g == n:
            steps = 0
            while g == 1 and steps < 1000:
                ys = (ys * ys + c) % n
                g = math.gcd(abs(x - ys), n)
                steps += 1
        if g != n:
            return g

def _factor_large(n, f):
    """Recursively fully factor composite/prime n into primes, into dict f."""
    if n <= 1:
        return
    if _is_prime_fast(n):
        f[n] = f.get(n, 0) + 1
        return
    d = _brent_rho(n)
    _factor_large(d, f)
    _factor_large(n // d, f)

def _factorize(n):
    """Full prime factorisation of n as a {prime: exponent} dict."""
    f = {}
    for q in _TRIAL:
        if q * q > n: break
        while n % q == 0:
            f[q] = f.get(q, 0) + 1
            n //= q
    if n > 1:
        _factor_large(n, f)
    return f

def _divs_of_sq(n):
    """All divisors of n^2, sorted ascending."""
    f = _factorize(n)
    ds = [1]
    for q, e in f.items():
        new, qe = [], 1
        for _ in range(2 * e):
            qe *= q
            for d in ds: new.append(d * qe)
        ds.extend(new)
    return sorted(ds)

def classify(p):
    """For p = 1 (mod 24): returns (code, A, d).
    code 0 = Case A, 1 = Case B NQR7, 2 = Case B QR7 solved, 3 = OPEN.
    Identical decision logic to session20_corrected_10B.py -- only the
    underlying factorisation calls changed."""
    m = (p + 3) // 4
    if any(q % 3 != 1 for q in _factorize(m)):
        return (0, 0, 0)                  # Case A: some factor = 2 (mod 3)
    if p % 7 in (3, 5, 6):
        return (1, 0, 0)                  # Case B, NQR mod 7: Prop nqr7
    for A in _TS:                         # Case B, QR mod 7: gateway search
        N = (p + A) >> 2
        B = p * N
        for d in _divs_of_sq(N):
            if d <= 1: continue
            if (B + d) % A: continue
            y = (B + d) // A
            if (B * y) % d: continue
            z = (B * y) // d
            if z > 0 and 4 * N * y * z == p * (y * z + N * z + N * y):
                return (2, A, d)
    return (3, 0, 0)

# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    import numpy as np
    try:
        import cupy as cp
        gpu = True
    except ImportError:
        cp = None
        gpu = False

    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    print("ERDOS-STRAUS -- SESSION 22: POLLARD-BRENT ENGINE")
    print("=" * 60)
    print(f"LIMIT = {LIMIT:,}   GPU sieve: {gpu}   workers: {WORKERS}")

    t0 = time.time()
    SMALL_PRIMES = _small_primes(SQRT_LIMIT)
    print(f"Sieve primes to {SQRT_LIMIT:,}: {len(SMALL_PRIMES):,}  ({time.time()-t0:.1f}s)")

    def sieve_segment(lo, hi):
        size = hi - lo + 1
        xp = cp if gpu else np
        s = xp.ones(size, dtype=xp.uint8)
        if lo == 0: s[0] = 0
        if lo <= 1 <= hi: s[1 - lo] = 0
        for q in SMALL_PRIMES:
            if q * q > hi: break
            start = max(q * q, ((lo + q - 1) // q) * q)
            if start > hi: continue
            s[start - lo::q] = 0
        return cp.asnumpy(s) if gpu else s

    def load_ckpt():
        if os.path.exists(CKPT):
            with open(CKPT) as f: return json.load(f)
        return None

    def save_ckpt(state):
        tmp = CKPT + ".tmp"
        with open(tmp, 'w') as f: json.dump(state, f)
        os.replace(tmp, CKPT)

    ck = load_ckpt()
    if ck:
        print(f"\nRESUMING from {ck['next_lo']:,}")
        start_lo   = ck['next_lo']
        totals     = ck['totals']
        Au         = Counter({int(k): v for k, v in ck['Au'].items()})
        mxA        = ck['mxA']
        opens      = ck['opens']
        witnesses  = ck['witnesses']
        milestones = ck['milestones']
    else:
        start_lo   = 0
        totals     = {'primes': 0, 'p2': 0, 'thm1': 0, 'thm2': 0, 'thm4': 0,
                      'caseA': 0, 'nqr7': 0, 'gw': 0}
        Au         = Counter()
        mxA        = 0
        opens      = []
        witnesses  = []
        milestones = {}

    pool = Pool(WORKERS, initializer=_init_worker)
    print(f"\nScanning {start_lo:,} -> {LIMIT:,} in segments of {SEG:,}")
    print("=" * 60)
    t_all = time.time()

    lo = start_lo
    while lo <= LIMIT:
        hi = min(lo + SEG - 1, LIMIT)
        t_seg = time.time()

        s = sieve_segment(lo, hi)
        primes = (np.nonzero(s)[0] + lo).astype(np.int64)
        n_seg = len(primes)

        r4  = primes % 4
        r8  = primes % 8
        r24 = primes % 24
        m_p2   = primes == 2
        m_thm1 = r4 == 3
        m_thm2 = (~m_thm1) & (r8 == 5)
        m_thm4 = (~m_thm1) & (~m_thm2) & (r24 == 17)
        m_m24  = (~m_p2) & (~m_thm1) & (~m_thm2) & (~m_thm4) & (r24 == 1)

        pend = [t for t in MILESTONES if lo <= t <= hi and str(t) not in milestones]
        base = {k: totals[k] for k in ('p2', 'thm1', 'thm2', 'thm4')}

        m24_list = [int(p) for p in primes[m_m24]]
        results = pool.map(classify, m24_list, chunksize=256)

        seg_gw = seg_open = 0
        it = iter(zip(m24_list, results))
        cur = next(it, None)
        run = {'caseA': totals['caseA'], 'nqr7': totals['nqr7'],
               'gw': totals['gw'], 'open': len(opens)}
        for t in sorted(pend) + [None]:
            while cur is not None and (t is None or cur[0] <= t):
                p, (code, A, d) = cur
                if code == 0:  run['caseA'] += 1
                elif code == 1: run['nqr7'] += 1
                elif code == 2:
                    run['gw'] += 1; seg_gw += 1
                    Au[A] += 1
                    if A > mxA: mxA = A
                    if A >= 100: witnesses.append((p, A, d))
                else:
                    run['open'] += 1; seg_open += 1
                    opens.append(p)
                cur = next(it, None)
            if t is not None:
                below = primes <= t
                milestones[str(t)] = {
                    'primes': int(totals['primes'] + int(below.sum())),
                    'p2':   int(base['p2']   + int((m_p2   & below).sum())),
                    'thm1': int(base['thm1'] + int((m_thm1 & below).sum())),
                    'thm2': int(base['thm2'] + int((m_thm2 & below).sum())),
                    'thm4': int(base['thm4'] + int((m_thm4 & below).sum())),
                    'caseA': run['caseA'], 'nqr7': run['nqr7'],
                    'gw': run['gw'], 'open': run['open'],
                    'mxA': mxA, 'Au': dict(Au),
                    'witnesses_geq100': [w for w in witnesses],
                }
                print(f"\n--- MILESTONE {t:,} ---")
                ms = milestones[str(t)]
                print(f"  primes={ms['primes']:,}  caseA={ms['caseA']:,}  "
                      f"nqr7={ms['nqr7']:,}  CaseB-QR7={ms['gw']:,}  "
                      f"open={ms['open']}  maxA={ms['mxA']}")
                print(f"  A-dist: {dict(sorted(Au.items()))}")
                if witnesses:
                    print(f"  witnesses A>=100: {witnesses}")
                print()

        totals['primes'] += n_seg
        totals['p2']   += int(m_p2.sum())
        totals['thm1'] += int(m_thm1.sum())
        totals['thm2'] += int(m_thm2.sum())
        totals['thm4'] += int(m_thm4.sum())
        totals['caseA'] = run['caseA']
        totals['nqr7']  = run['nqr7']
        totals['gw']    = run['gw']

        print(f"[{lo//SEG+1:3d}] {lo:>14,}-{hi:>14,}  {n_seg:,} primes  "
              f"gw={seg_gw}  open={seg_open}  maxA={mxA}  "
              f"{time.time()-t_seg:.0f}s")

        save_ckpt({'next_lo': hi + 1, 'totals': totals, 'Au': dict(Au),
                   'mxA': mxA, 'opens': opens[:1000],
                   'witnesses': witnesses, 'milestones': milestones})
        lo = hi + 1

    pool.close(); pool.join()

    print("\n" + "=" * 60)
    print("SCAN COMPLETE")
    print("=" * 60)
    tot = totals['primes']
    print(f"  Wall: {(time.time()-t_all)/3600:.2f}h   Primes: {tot:,}")
    labels = [("p = 2", 'p2'), ("p = 3 (mod 4)  [A=1]", 'thm1'),
              ("p = 5 (mod 8)  [A=3,d=2]", 'thm2'),
              ("p = 17 (mod 24) [A=3,d=p]", 'thm4'),
              ("Case A  [A=3,d=q]", 'caseA'),
              ("Case B NQR7  [A=7]", 'nqr7'),
              ("Case B QR7  [gateway]", 'gw')]
    for lab, k in labels:
        print(f"  {lab:<28} {totals[k]:>13,}  {100.0*totals[k]/tot:7.3f}%")
    print(f"  {'OPEN':<28} {len(opens):>13}")
    print(f"\n  max A = {mxA}")
    print(f"  A-dist: {dict(sorted(Au.items()))}")
    print(f"  witnesses A>=100: {witnesses}")

if __name__ == '__main__':
    main()
