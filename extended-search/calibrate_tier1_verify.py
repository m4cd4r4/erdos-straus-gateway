"""Calibration: cost of exhaustively re-verifying every A>=100 witness in the
extended-search range, independently of the primary engine.

For each sampled witness (p, A, d):
  1. Validate the witness itself: p prime, d | N_A^2, d = t_A mod A.
  2. Exhaustively check every candidate gateway B = 3 mod 4, B < A: compute
     N_B = (p+B)/4, factorise it, build the full divisor-residue set mod B,
     and confirm t_B is NOT in it (i.e. B genuinely fails). This is the
     referee's method for the two original records, applied to a random
     sample instead of just those two.

Times the sample, extrapolates to the full 19,540 A>=100 witnesses, and
stratifies by A so the extrapolation isn't misled by an unrepresentative
sample of a distribution that is skewed toward small A (median 107).

Run: python calibrate_tier1_verify.py [N_SAMPLE]
"""
import sys, json, time, random
from sympy import factorint, primerange, isprime

CKPT = 'results/session23_checkpoint.json'
N_SAMPLE = int(sys.argv[1]) if len(sys.argv) > 1 else 100

GW_CACHE = {}


def gateways_below(A):
    if A not in GW_CACHE:
        GW_CACHE[A] = [q for q in primerange(3, A) if q % 4 == 3]
    return GW_CACHE[A]


def divisor_residues(facts, B):
    D = {1 % B}
    for q, e in facts.items():
        if q % B == 0:
            continue
        qm, pw, new = q % B, 1, set()
        for j in range(2 * e + 1):
            if j:
                pw = (pw * qm) % B
            for r in D:
                new.add((r * pw) % B)
        D = new
    return D


def verify_witness(p, A, d):
    """Full exhaustive verification of one witness. Returns (ok, n_gateways_checked)."""
    assert isprime(p), p
    N_A = (p + A) // 4
    assert (N_A * N_A) % d == 0, (p, A, d, 'd does not divide N_A^2')
    t_A = (-4 * N_A * N_A) % A
    assert d % A == t_A, (p, A, d, 'd != t_A mod A')

    n_checked = 0
    for B in gateways_below(A):
        n_checked += 1
        N_B = (p + B) // 4
        if N_B % B == 0:
            continue
        t_B = (-4 * N_B * N_B) % B
        facts = factorint(N_B)
        D = divisor_residues(facts, B)
        assert t_B not in D, (p, A, B, 'SMALLER GATEWAY ALSO SUCCEEDS - claimed A is not least!')
    return True, n_checked


def main():
    d = json.load(open(CKPT))
    witnesses = [w for w in d['witnesses'] if w[1] >= 100]
    print(f'population: {len(witnesses)} witnesses with A >= 100')

    random.seed(20260728)
    sample = random.sample(witnesses, min(N_SAMPLE, len(witnesses)))
    sample.sort(key=lambda w: w[1])

    print(f'sampling {len(sample)} at random, verifying each exhaustively...')
    print(f"  {'p':>14} {'A':>5} {'gateways checked':>16} {'time (s)':>9}")

    results = []
    t_start = time.time()
    for p, A, dd in sample:
        t0 = time.time()
        ok, n = verify_witness(p, A, dd)
        dt = time.time() - t0
        results.append((A, n, dt))
        print(f'  {p:>14} {A:>5} {n:>16} {dt:>9.3f}')
    total_time = time.time() - t_start

    print()
    print(f'sample total: {total_time:.1f}s for {len(sample)} witnesses '
          f'({total_time/len(sample):.3f}s/witness average)')

    # Stratify by A-bucket, since cost scales with gateway count ~ A.
    buckets = [(100, 150), (150, 250), (250, 400), (400, 500)]
    print()
    print('cost by A-bucket (in the sample):')
    print('  ' + f"{'A range':>12} {'n in sample':>12} {'mean time (s)':>14}")
    bucket_stats = {}
    for lo, hi in buckets:
        in_b = [t for A, n, t in results if lo <= A < hi]
        if in_b:
            bucket_stats[(lo, hi)] = sum(in_b) / len(in_b)
            print('  ' + f'{f"{lo}-{hi}":>12} {len(in_b):>12} {sum(in_b)/len(in_b):>14.3f}')

    # Extrapolate using the population's real A-distribution, not the flat
    # sample average, so the projection isn't skewed by how the sample landed.
    pop_bucket_counts = {b: 0 for b in buckets}
    for w in witnesses:
        for lo, hi in buckets:
            if lo <= w[1] < hi:
                pop_bucket_counts[(lo, hi)] += 1
                break

    print()
    print('extrapolation to the full population, using per-bucket sample cost:')
    print('  ' + f"{'A range':>12} {'population n':>13} {'proj. core-seconds':>19}")
    total_proj = 0.0
    for b in buckets:
        if b in bucket_stats:
            proj = pop_bucket_counts[b] * bucket_stats[b]
            total_proj += proj
            print('  ' + f'{f"{b[0]}-{b[1]}":>12} {pop_bucket_counts[b]:>13} {proj:>19.1f}')
        elif pop_bucket_counts[b]:
            print('  ' + f'{f"{b[0]}-{b[1]}":>12} {pop_bucket_counts[b]:>13} '
                  f'{"NO SAMPLE - projection unreliable":>19}')

    print()
    print(f'TOTAL projected: {total_proj:.0f} core-seconds = {total_proj/3600:.2f} core-hours')
    for workers in (10, 16):
        print(f'  wall time at {workers} workers: {total_proj/3600/workers:.2f} h')

    json.dump(dict(sample=[(p, A, dd) for p, A, dd in sample], results=results,
                   total_projected_core_seconds=total_proj),
              open('results/calibrate_tier1_results.json', 'w'))


if __name__ == '__main__':
    main()
