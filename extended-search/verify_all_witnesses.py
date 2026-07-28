"""Tier 1 verification, full population: every A>=100 witness, exhaustively.

Same method as calibrate_tier1_verify.py's verify_witness(), applied to the
entire witness list rather than a random sample. See tier1_verification.md
for the result (19,540 / 19,540 pass, 55.8s, max A confirmed at 479) and why
this was cheap enough to just run rather than keep calibrating.

Run: python verify_all_witnesses.py
"""
import json, time
from sympy import factorint, primerange, isprime

CKPT = 'results/session23_checkpoint.json'
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
    assert isprime(p), p
    N_A = (p + A) // 4
    assert (N_A * N_A) % d == 0, (p, A, d, 'd does not divide N_A^2')
    t_A = (-4 * N_A * N_A) % A
    assert d % A == t_A, (p, A, d, 'd != t_A mod A')
    for B in gateways_below(A):
        N_B = (p + B) // 4
        if N_B % B == 0:
            continue
        t_B = (-4 * N_B * N_B) % B
        if t_B in divisor_residues(factorint(N_B), B):
            return False, B
    return True, None


def main():
    d = json.load(open(CKPT))
    witnesses = [w for w in d['witnesses'] if w[1] >= 100]
    print(f'verifying ALL {len(witnesses)} witnesses A >= 100, exhaustively...')

    t0 = time.time()
    n_ok, fails = 0, []
    for p, A, dd in witnesses:
        ok, culprit = verify_witness(p, A, dd)
        if ok:
            n_ok += 1
        else:
            fails.append((p, A, dd, culprit))
    elapsed = time.time() - t0

    print(f'done in {elapsed:.1f}s')
    print(f'OK: {n_ok} / {len(witnesses)}')
    print(f'FAILURES: {len(fails)}')
    for f in fails[:20]:
        print('  ', f)
    print(f'max A independently confirmed: {max(A for _, A, _ in witnesses)}')

    json.dump(dict(n_witnesses=len(witnesses), n_ok=n_ok, failures=fails,
                    elapsed_s=elapsed,
                    max_A_confirmed=max(A for _, A, _ in witnesses)),
              open('results/tier1_full_verification.json', 'w'))


if __name__ == '__main__':
    main()
