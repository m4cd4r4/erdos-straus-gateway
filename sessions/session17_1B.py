#!/usr/bin/env python3
"""
ERDOS-STRAUS -- SESSION 17: SCALE TO 1,000,000,000
====================================================
Session 16: 100% at 100M. Max A needed: 239 (1 prime escaped A<200).
Session 17: Push to 1B. Does max A stay bounded?
Memory-optimized: trial division, streaming scan, free sieve after use.
"""
import sys, time
from math import isqrt
from collections import Counter

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
print("ERDOS-STRAUS -- SESSION 17: SCALE TO 1,000,000,000")
print("=" * 60)

LIMIT = 1_000_000_000

t0 = time.time()
print("Sieving primes to 1B (~1GB RAM)...")
ip = bytearray([1]) * (LIMIT + 1)
ip[0] = ip[1] = 0
for i in range(2, isqrt(LIMIT) + 1):
    if ip[i]:
        ip[i*i::i] = bytearray(len(ip[i*i::i]))
print("Sieve complete (%.1fs)" % (time.time() - t0))

SP = [i for i in range(2, 100000) if ip[i]]
PS = set(SP)

def factorize(n):
    if n <= 1: return {}
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

def chk(p, t, d):
    A = 3 + 4*t
    xn = p + A
    if xn % 4: return False
    x = xn // 4
    B = p * x
    if (B + d) % A: return False
    y = (B + d) // A
    By = B * y
    if By % d: return False
    z = By // d
    return z > 0 and 4*x*y*z == p*(y*z + x*z + x*y)

TS0 = [t for t in range(1, 60) if (3+4*t) in PS and (3+4*t) < 200]
TS1 = [t for t in range(1, 260) if (3+4*t) in PS and 200 <= (3+4*t) < 1000]
TS2 = [t for t in range(1, 2600) if (3+4*t) in PS and 1000 <= (3+4*t) < 10000]
print("A values: %d (<200), +%d (<1K), +%d (<10K)" % (len(TS0), len(TS1), len(TS2)))

print("\n" + "=" * 60)
print("STREAMING SCAN")
print("=" * 60)

t0 = time.time()
c = [0]*7  # p2, thm1, thm2, thm4, thm5, thm7, gw
n_ext = 0; opens = []; Au = Counter(); dgt = 0; mxA = 0; tot = 0; cbq = 0

for p in range(2, LIMIT + 1):
    if not ip[p]: continue
    tot += 1
    if p == 2: c[0] = 1; continue
    r4 = p % 4
    if r4 == 3: c[1] += 1; continue
    if p % 8 == 5: c[2] += 1; continue
    if p % 24 == 17: c[3] += 1; continue
    if p % 24 != 1: continue
    tmp = (p+3)//4; ok = True; d2 = 2; n2 = tmp
    while d2*d2 <= n2:
        if n2 % d2 == 0:
            if d2 % 3 != 1: ok = False; break
            while n2 % d2 == 0: n2 //= d2
        d2 += 1
    if ok and n2 > 1 and n2 % 3 != 1: ok = False
    if not ok: c[4] += 1; continue
    if p % 7 in (3, 5, 6): c[5] += 1; continue
    cbq += 1; found = False
    for ts_list, phase in [(TS0, 0), (TS1, 1), (TS2, 2)]:
        for t in ts_list:
            A = 3+4*t
            if (p+A) % 4: continue
            Nt = (p+A)//4
            for d in divs_of_sq(Nt):
                if d <= 1: continue
                if chk(p, t, d):
                    found = True; Au[A] += 1
                    if A > mxA: mxA = A
                    if d > Nt: dgt += 1
                    if phase == 0: c[6] += 1
                    else: n_ext += 1
                    break
            if found: break
        if found: break
    if not found: opens.append(p)
    if cbq % 200000 == 0:
        el = time.time()-t0
        ng = c[6]+n_ext
        print("  %s CaseB-QR7 | %s solved | %d open | maxA=%d | %.0fs" % (
            format(cbq, ","), format(ng, ","), len(opens), mxA, el))

del ip
elapsed = time.time()-t0
ngw = c[6]+n_ext

print("\nScan complete (%.0fs)" % elapsed)
print("  Total primes: %s" % format(tot, ","))
print("  Case B QR7: %s" % format(cbq, ","))
print("  Solved (A<200): %s" % format(c[6], ","))
print("  Solved (A>=200): %s" % format(n_ext, ","))
print("  Open: %d" % len(opens))
print("  Max A: %d" % mxA)

print("\n" + "=" * 60)
print("A-VALUE DISTRIBUTION")
print("=" * 60)
print("\n  %6s  %10s  %8s  %8s" % ("A", "count", "%", "cumul%"))
cum = 0
for A, cnt in Au.most_common(40):
    cum += cnt
    print("  %6d  %10s  %7.2f%%  %7.2f%%" % (A, format(cnt,","), 100.0*cnt/max(ngw,1), 100.0*cum/max(ngw,1)))
print("\nd > N_t (N_t^2 needed): %s / %s = %.2f%%" % (format(dgt,","), format(ngw,","), 100.0*dgt/max(ngw,1)))

if opens:
    print("\n%d OPEN primes:" % len(opens))
    print("First 20:", opens[:20])

print("\n" + "=" * 60)
print("COVERAGE TABLE (PRIMES TO 1B)")
print("=" * 60)
labels = ["p=2","p=3(mod 4) [Thm 1]","p=5(mod 8) [Thm 2]","p=17(mod 24) [Thm 4]",
          "Case A [Thm 5]","Case B NQR-mod-7 [Thm 7]","Case B QR-mod-7 [Gateway]"]
vals = [c[0],c[1],c[2],c[3],c[4],c[5],ngw]
proven = sum(vals)
print()
for lab, v in zip(labels, vals):
    print("  %-48s  %12s  %7.3f%%" % (lab, format(v,","), 100.0*v/tot))
print("  " + "-"*48 + "  " + "-"*12 + "  " + "-"*8)
print("  %-48s  %12s  %7.4f%%" % ("PROVEN TOTAL", format(proven,","), 100.0*proven/tot))
print("  %-48s  %12s  %7.5f%%" % ("OPEN", str(len(opens)), 100.0*len(opens)/tot))

print("\nSCALE PROGRESSION:")
print("  10^6:  78,498 primes       100%% (maxA=199)")
print("  10^7:  664,579 primes      100%% (maxA=199)")
print("  10^8:  5,761,455 primes    100%% (maxA=239)")
print("  10^9:  %s primes  %.4f%% (maxA=%d)" % (format(tot,","), 100.0*proven/tot, mxA))

if not opens:
    print("\n" + "=" * 60)
    print("ALL PRIMES TO 1,000,000,000 PROVEN ALGEBRAICALLY!")
    print("=" * 60)
