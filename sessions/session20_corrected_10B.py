#!/usr/bin/env python3
"""
ERDOS-STRAUS -- SESSION 20: CORRECTED RE-VERIFICATION TO 10^10
==============================================================
Re-runs the 10^10 verification with the CORRECT Case A / Case B split.

Why this session exists (2026-07-05 referee finding):
  session18_10B.py's `is_case_a()` returned True when ALL prime factors of
  (p+3)/4 are congruent to 1 (mod 3) -- but that is Case B, the hard class.
  Session 18 therefore counted every hard prime as "proven by Thm 5" and
  skipped it, while gateway-searching the complementary easy class
  (Case A intersect QR7, which Prop mod3(c) already covers with A = 3).
  Sessions 15-17 and paper/collect_stats.py used the correct split, so
  results to 10^9 stand (maxA = 239 at 10^8 and 10^9 per sessions 16/17).
  Session 18's 10^10 dataset (18,189,229 "hard" primes, maxA 251, the
  25-value A dictionary) describes the wrong set and is void.

Classification (correct, matches sessions 15-17 and the paper):
  p = 2                                -> p2
  p ≡ 3 (mod 4)                        -> thm1  (classical, A = 1)
  p ≡ 5 (mod 8)                        -> thm2  (A = 3, d = 2)
  p ≡ 17 (mod 24)                      -> thm4  (A = 3, d = p)
  p ≡ 1 (mod 24), some q|N, q ≢ 1 (3)  -> caseA (A = 3, d = q; Prop mod3(c))
  p ≡ 1 (mod 24), Case B, NQR mod 7    -> nqr7  (A = 7; Prop nqr7)
  p ≡ 1 (mod 24), Case B, QR mod 7     -> gateway search over prime A ≡ 3 (mod 4)

Usage:
  python session20_corrected_10B.py [LIMIT] [CKPT_PATH]
  Defaults: LIMIT = 10^10, CKPT = results/session20_checkpoint.json
  Resumable: checkpoints after every segment.

Baselines the run must reproduce (correct-classifier ground truth):
  10^6: gw = 2,269, maxA = 79   (paper/stats_1M.json, collect_stats.py)
  10^8: gw = 146,251, maxA = 239, one prime with A >= 200  (session 16)
  10^9: gw = 1,212,383, maxA = 239                         (session 17)
"""
import sys, time, json, os
from math import isqrt
from collections import Counter
from multiprocessing import Pool

LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000_000_000
CKPT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "results", "session20_checkpoint.json")

SEG        = 100_000_000
SQRT_LIMIT = max(isqrt(LIMIT), 1000)
WORKERS    = 16
MILESTONES = [10**6, 10**7, 10**8, 10**9, 10**10, 10**11]

# ── Worker globals (built once per process via initializer) ─────────────────

_SP = None      # primes to sqrt(LIMIT), for factorisation
_TS = None      # candidate A values: primes ≡ 3 (mod 4), ascending, < 10^4

def _small_primes(n):
    ip = bytearray([1]) * (n + 1)
    ip[0] = ip[1] = 0
    for i in range(2, isqrt(n) + 1):
        if ip[i]:
            ip[i*i::i] = bytearray(len(ip[i*i::i]))
    return [i for i in range(2, n + 1) if ip[i]]

def _init_worker(sqrt_limit):
    global _SP, _TS
    _SP = _small_primes(sqrt_limit)
    _TS = [A for A in _small_primes(10_000) if A % 4 == 3]

def _divs_of_sq(n):
    """All divisors of n^2, from the factorisation of n.
    Valid while n < sqrt_limit^2: the cofactor after trial division is prime."""
    f = {}
    m = n
    for q in _SP:
        if q * q > m: break
        while m % q == 0:
            f[q] = f.get(q, 0) + 1
            m //= q
    if m > 1: f[m] = f.get(m, 0) + 1
    ds = [1]
    for q, e in f.items():
        new, qe = [], 1
        for _ in range(2 * e):
            qe *= q
            for d in ds: new.append(d * qe)
        ds.extend(new)
    return sorted(ds)

def classify(p):
    """For p ≡ 1 (mod 24): returns (code, A, d).
    code 0 = Case A, 1 = Case B NQR7, 2 = Case B QR7 solved, 3 = OPEN."""
    m = (p + 3) // 4
    for q in _SP:
        if q * q > m: break
        if m % q == 0:
            if q % 3 != 1:
                return (0, 0, 0)          # Case A: q ≡ 2 (mod 3) divides N
            while m % q == 0: m //= q
    if m > 1 and m % 3 != 1:
        return (0, 0, 0)                  # large prime cofactor ≡ 2 (mod 3)
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
    print("ERDOS-STRAUS -- SESSION 20: CORRECTED RE-VERIFICATION")
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
        witnesses  = []          # (p, A, d) for every solution with A >= 100
        milestones = {}

    pool = Pool(WORKERS, initializer=_init_worker, initargs=(SQRT_LIMIT,))
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

        print(f"[{lo//SEG+1:3d}] {lo:>13,}-{hi:>13,}  {n_seg:,} primes  "
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
    labels = [("p = 2", 'p2'), ("p ≡ 3 (mod 4)  [A=1]", 'thm1'),
              ("p ≡ 5 (mod 8)  [A=3,d=2]", 'thm2'),
              ("p ≡ 17 (mod 24) [A=3,d=p]", 'thm4'),
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
