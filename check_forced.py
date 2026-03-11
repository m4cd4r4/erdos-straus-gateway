#!/usr/bin/env python3
"""Check which candidate A values have forced prime divisors of N_A."""
from math import isqrt

def is_prime(n):
    if n < 2: return False
    for d in range(2, isqrt(n)+1):
        if n % d == 0: return False
    return True

candidates = [a for a in range(3, 240) if a % 4 == 3 and is_prime(a)]
print(f"{len(candidates)} candidates: {candidates[:5]}...{candidates[-3:]}")

QR_mod = {}
for A in candidates:
    QR_mod[A] = {x*x % A for x in range(1, A)} - {0}

def forced_prime_factors(A):
    """Find primes q that always divide N_A = (p+A)/4 for Case B QR7 primes p.

    p ≡ 1 mod 24. N_A = (p+A)/4. For q to always divide N_A,
    we need q | (p+A)/4 for all p ≡ 1 mod 24.
    This means (1+A)/4 ≡ 0 mod q (i.e., q | (1+A)/4).
    """
    val = (1 + A) // 4  # N_A when p=1 (the residue mod 24)
    result = []
    for q in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
        if val % q == 0:
            # q | (1+A)/4. Does this mean q | N_A for ALL p ≡ 1 mod 24?
            # N_A = (p+A)/4 ≡ (1+A)/4 mod q. If q | (1+A)/4 then yes.
            result.append(q)
    return result

print("\nForced prime divisors of N_A for Case B QR7 primes (p = 1 mod 24):")
sharp_candidates = []
for A in candidates:
    ff = forced_prime_factors(A)
    if ff:
        # Check which forced factors are QR mod A and generate Q_A
        useful = []
        for g in ff:
            if g % A not in QR_mod[A]:
                continue  # g must be QR mod A
            # Check if g generates Q_A
            order = 1
            cur = g % A
            while cur != 1:
                cur = (cur * g) % A
                order += 1
            if order == len(QR_mod[A]):
                # g generates Q_A; check if v_g(N_A) sufficient
                # min v_g((1+A)/4) tells us the guaranteed valuation
                val = (1 + A) // 4
                v = 0
                while val % g == 0:
                    v += 1
                    val //= g
                useful.append((g, order, v, 2*v >= len(QR_mod[A])-1))
        if useful:
            print(f"  A={A:3d}: forced={ff}, useful={useful}")
            if any(u[3] for u in useful):
                sharp_candidates.append(A)
        else:
            print(f"  A={A:3d}: forced={ff} (none generate Q_A or insufficient valuation)")

print(f"\nPotentially sharp candidates: {sharp_candidates}")
