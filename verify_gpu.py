#!/usr/bin/env python3
"""
GPU-accelerated Erdos-Straus verification using CuPy RawKernel.

RTX 5000 (Turing, 3072 cores, 16GB VRAM).

Strategy:
  CPU  : segmented sieve (proven correct, same logic as verify_segmented.py)
  CPU  : algebraic filters (five mod checks eliminate ~97.6% of primes)
  GPU  : gateway kernel via CuPy RawKernel (one thread per Case B QR7 prime)

Expected speedup over single-core CPU: 100-500x on gateway checks.
Expected total runtime to 10^10: ~30-90 minutes vs 12-20 hours.

Usage:
  python verify_gpu.py 10000000000          # verify to 10^10
  python verify_gpu.py 1000000000           # quick test to 10^9 (~2 min)
  python verify_gpu.py 10000000000 512      # custom threads-per-block
"""

import sys, time
from math import isqrt

# Force line-buffered stdout so progress appears immediately in log files
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
import cupy as cp

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

LIMIT      = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000_000_000
TPB        = int(sys.argv[2]) if len(sys.argv) > 2 else 256
SEGMENT    = 500_000_000   # sieve segment (500MB RAM)
BATCH_SIZE = 2_000_000     # GPU batch size
MAX_A      = 239
SQRT_LIMIT = isqrt(LIMIT)

# -----------------------------------------------------------------------------
# Precompute small primes and candidate A values
# -----------------------------------------------------------------------------

def _sieve_small(n):
    ip = bytearray([1]) * (n + 1)
    ip[0] = ip[1] = 0
    for i in range(2, isqrt(n) + 1):
        if ip[i]:
            ip[i*i::i] = bytearray(len(ip[i*i::i]))
    return [i for i in range(2, n + 1) if ip[i]]

print(f"Sieving small primes to {SQRT_LIMIT:,}... ", end="", flush=True)
small_primes = _sieve_small(SQRT_LIMIT)
print(f"{len(small_primes):,} found")

candidate_as = [a for a in small_primes if a % 4 == 3 and a <= MAX_A]
cand_np      = np.array(candidate_as, dtype=np.int64)
N_CANDS      = len(candidate_as)
print(f"Candidate A values ({N_CANDS}): {candidate_as[:5]} ... {candidate_as[-3:]}")

# -----------------------------------------------------------------------------
# CuPy RawKernel  –  gateway check (one CUDA thread per prime)
# -----------------------------------------------------------------------------
# For each prime p the thread:
#   1. Loops over candidate A values in order.
#   2. For the first A where 4|(p+A), factorises N=(p+A)/4 by trial division.
#   3. Builds the set of residues of divisors of N^2 mod A iteratively.
#   4. If the target  t = (-4N^2) mod A  is in that set, writes A to results[].
# results[i] = smallest A that solved prime i, or -1 for an open case.

_KERNEL_SRC = r"""
extern "C" __global__
void gateway_kernel(
    const long long* __restrict__ primes,   // Case B QR7 primes
    const long long* __restrict__ cands,    // candidate A values
    int n_cands,
    long long* __restrict__ results,
    int n_primes)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_primes) return;

    long long p = primes[idx];
    results[idx] = -1LL;   // default: open case

    // Thread-local storage
    long long fac[12];
    int       exp_a[12];
    char      res[240];    // residue set, size >= MAX_A = 239
    char      nxt[240];

    for (int ai = 0; ai < n_cands; ++ai) {
        long long A = cands[ai];

        if ((p + A) % 4LL != 0LL) continue;

        long long N = (p + A) / 4LL;

        // -- Trial-divide N --------------------------------------
        int nf  = 0;
        long long rem = N;
        long long d   = 2LL;
        while (d * d <= rem && nf < 12) {
            if (rem % d == 0LL) {
                fac[nf]   = d;
                exp_a[nf] = 0;
                while (rem % d == 0LL) { ++exp_a[nf]; rem /= d; }
                ++nf;
            }
            ++d;
        }
        if (rem > 1LL && nf < 12) {
            fac[nf]   = rem;
            exp_a[nf] = 1;
            ++nf;
        }

        // -- Target residue: t = (-4N^2) mod A -------------------
        long long Nm     = N % A;
        long long target = (A - (4LL * Nm * Nm) % A) % A;

        // -- Initialise residue set = {1} ------------------------
        int Ai = (int)A;
        for (int r = 0; r < Ai; ++r) res[r] = 0;
        res[1] = 1;

        bool found = (target == 1LL);

        // -- Build divisors of N^2 mod A iteratively -------------
        for (int fi = 0; fi < nf && !found; ++fi) {
            for (int r = 0; r < Ai; ++r) nxt[r] = 0;

            long long q_pow = 1LL;
            int max_k = 2 * exp_a[fi];
            for (int k = 0; k <= max_k; ++k) {
                for (int r = 0; r < Ai; ++r) {
                    if (res[r]) {
                        long long nr = ((long long)r * q_pow) % A;
                        nxt[nr] = 1;
                        if (nr == target) found = true;
                    }
                }
                q_pow = (q_pow * fac[fi]) % A;
            }
            for (int r = 0; r < Ai; ++r) res[r] = nxt[r];
        }

        if (found) {
            results[idx] = A;
            return;
        }
    }
    // results[idx] remains -1  (open case)
}
"""

print("Compiling CUDA kernel... ", end="", flush=True)
_kernel = cp.RawKernel(_KERNEL_SRC, "gateway_kernel")
cand_gpu = cp.array(cand_np, dtype=cp.int64)
print("done")

# -----------------------------------------------------------------------------
# CPU algebraic filters  (identical to verify_segmented.py)
# -----------------------------------------------------------------------------

def _prop_3mod4(p):   return p % 4 == 3
def _prop_5mod8(p):   return p % 8 == 5
def _prop_17mod24(p): return p % 24 == 17

def _prop_nqr7(p):
    return p % 24 == 1 and p % 7 in (3, 5, 6)

def _prop_caseA(p):
    if p % 24 != 1:
        return False
    n = (p + 3) // 4
    d = 2
    while d * d <= n:
        while n % d == 0:
            if d % 3 == 2:
                return True
            n //= d
        d += 1
    return n > 1 and n % 3 == 2

# -----------------------------------------------------------------------------
# GPU batch runner
# -----------------------------------------------------------------------------

def _flush_batch(batch, stats):
    if not batch:
        return
    n        = len(batch)
    p_gpu    = cp.array(batch, dtype=cp.int64)
    res_gpu  = cp.full(n, -1, dtype=cp.int64)
    blocks   = (n + TPB - 1) // TPB
    _kernel(
        (blocks,), (TPB,),
        (p_gpu, cand_gpu, np.int32(N_CANDS), res_gpu, np.int32(n))
    )
    results = res_gpu.get()
    for i, r in enumerate(results):
        if r == -1:
            stats['open'] += 1
            print(f"\n  *** OPEN: p={batch[i]} not solved by any A <= {MAX_A} ***",
                  flush=True)
        elif int(r) > stats['max_A']:
            stats['max_A']   = int(r)
            stats['hardest'] = int(batch[i])

# -----------------------------------------------------------------------------
# Warmup / correctness check
# -----------------------------------------------------------------------------

def _warmup():
    """Verify known results using genuine Case B QR7 primes."""
    # These pass all five algebraic filters; 32349601 requires A=239.
    tests  = [193, 673, 2521, 32349601]
    batch  = np.array(tests, dtype=np.int64)
    p_gpu  = cp.array(batch, dtype=cp.int64)
    r_gpu  = cp.full(len(tests), -1, dtype=cp.int64)
    _kernel((1,), (32,), (p_gpu, cand_gpu, np.int32(N_CANDS), r_gpu, np.int32(len(tests))))
    r      = r_gpu.get()
    print(f"Warmup check (primes {tests}):")
    ok = True
    for p, a in zip(tests, r):
        print(f"  p={p:>12,}  A={a}")
        if a == -1:
            ok = False
            print("  *** FAIL: expected a solution ***")
    return ok

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    t0  = time.time()
    dev = cp.cuda.runtime.getDeviceProperties(0)
    print(f"\nErdos-Straus GPU Verification   limit = {LIMIT:,}")
    print(f"Device  : {dev['name'].decode()}")
    print(f"TPB     : {TPB}   Segment : {SEGMENT:,}   Batch : {BATCH_SIZE:,}")
    print("=" * 70)

    if not _warmup():
        sys.exit("Warmup FAILED - check kernel correctness before long run")
    print()

    stats = {'max_A': 7, 'hardest': None, 'open': 0}
    total_primes = 0
    total_caseB  = 0
    pending      = []
    seg_start    = 2
    seg_num      = 0

    while seg_start <= LIMIT:
        seg_end  = min(seg_start + SEGMENT - 1, LIMIT)
        seg_size = seg_end - seg_start + 1

        # -- Segmented CPU sieve ---------------------------------------------
        sieve = bytearray([1]) * seg_size
        for sp in small_primes:
            start_idx = ((seg_start + sp - 1) // sp) * sp - seg_start
            if seg_start <= sp <= seg_end:
                sieve[sp - seg_start] = 1
                start_idx = sp * 2 - seg_start
            if 0 <= start_idx < seg_size:
                sieve[start_idx::sp] = bytearray(len(sieve[start_idx::sp]))

        # -- Algebraic filters + collect Case B QR7 primes ------------------
        for i in range(seg_size):
            p = seg_start + i
            if not sieve[i] or p < 3:
                continue
            total_primes += 1
            if (_prop_3mod4(p) or _prop_5mod8(p) or _prop_17mod24(p)
                    or _prop_caseA(p) or _prop_nqr7(p)):
                continue
            total_caseB += 1
            pending.append(p)

        # -- Flush GPU batch -------------------------------------------------
        if len(pending) >= BATCH_SIZE or seg_end >= LIMIT:
            _flush_batch(pending, stats)
            pending = []

        seg_num += 1
        elapsed  = time.time() - t0
        pct      = 100.0 * seg_end / LIMIT
        rate     = seg_end / elapsed if elapsed > 0 else 1
        eta      = (LIMIT - seg_end) / rate

        if seg_num % 5 == 0 or seg_end >= LIMIT:
            print(f"  [{pct:5.1f}%] to {seg_end:13,} | "
                  f"CaseB: {total_caseB:9,} | maxA: {stats['max_A']:3d} | "
                  f"{elapsed:6.0f}s elapsed | ETA {eta/3600:.2f}h",
                  flush=True)

        seg_start = seg_end + 1

    elapsed = time.time() - t0
    print()
    print("=" * 70)
    print(f"RESULT: {total_primes:,} primes verified up to {LIMIT:,}")
    print(f"  Case B QR7 primes   : {total_caseB:,}")
    print(f"  Maximum A required  : {stats['max_A']}")
    print(f"  Hardest prime found : {stats['hardest']}")
    print(f"  Open cases          : {stats['open']}")
    print(f"  Total time          : {elapsed/3600:.2f} hours")

    if stats['open'] == 0:
        tag = "VERIFIED" if stats['max_A'] <= MAX_A else "*** NEW MAX A ***"
        print(f"\n{tag}: Erdos-Straus holds for all primes p <= {LIMIT:,}")
        if stats['hardest']:
            print(f"  New hardest: p={stats['hardest']}, requires A={stats['max_A']}")
            print(f"  Prev hardest (to 10^9): p=32,349,601, A=239")

if __name__ == "__main__":
    main()
