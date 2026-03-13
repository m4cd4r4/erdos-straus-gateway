# Erdos-Straus Conjecture: Paper Plan (Revised March 2026)

## Working Title

**"Bounded Gateway Parameters for the Erdos-Straus Conjecture"**

*Subtitle: A Structural Analysis of Egyptian Fraction Decompositions to 10^10*

---

## Core Thesis

The paper's headline is NOT the computational verification (10^18 already done
by Mihnea-Dumitru 2025). The headline is:

> **The bounded-A conjecture**: the maximum gateway parameter A needed to solve
> any prime p <= N grows extremely slowly - by only ~52 across five decades of
> prime range - and may be absolutely bounded.

This structural observation appears new in the literature.

---

## Key Simplifications from Audit

### 1. ONE theorem, not fourteen
Sessions 7-12 evolved through single/semiprime/triple gateways, but these are all
instances of ONE parametric identity. The paper states it once.

### 2. The parameter t is eliminated
Replace t with A directly. There is no mathematical content in t = (A-3)/4.

### 3. The final identity check is a tautology
The verification 4xyz = p(yz + xz + xy) holds AUTOMATICALLY whenever the
integrality conditions are satisfied. Proved algebraically (see S3 proof).

### 4. Theorems 2, 4, and 5 merge into one
All three exploit "B = pN has a factor q not-equiv 1 (mod 3)" at A = 3.

### 5. The 7-layer hierarchy reduces to 4 layers
  1. p equiv 3 (mod 4): classical identity (~50%)
  2. p equiv 1 (mod 4), pN has factor q not-equiv 1 (mod 3) at A=3: unified (~47.5%)
  3. p equiv 1 (mod 24), Case B, p NQR mod 7: gateway at A=7 (~1.25%)
  4. p equiv 1 (mod 24), Case B, p QR mod 7: gateway search A < 251 (~1.25%)

### 6. The residue condition has a clean single statement
Instead of two checks (A | B+d) and (d | By), the condition is:
  d equiv -p^2 * 4^{-1} (mod A)  with  d | N^2

---

## Positioning & Contribution

### What's known (updated March 2026)
- Conjecture verified computationally to **10^18** (Mihnea & Dumitru, 2025;
  arXiv:2509.00128), using Salez's modular filtering algorithm
- Prior records: 10^17 (Salez, 2014), 10^14 (Swett, ~2000)
- Elsholtz-Tao (2013): solution count f(p) ~ (log p)^3 on average;
  Type I/II solution classification
- Tao: holds for "almost all" n; QR obstruction prevents covering systems from
  eliminating the Mordell subset
- Mordell subset: p mod 840 in {1, 121, 169, 289, 361, 529} are the "hard" primes
- 804/840 residue classes have explicit formulas; 36 remain
- arXiv:2403.16047 reformulates the hard case as (x, d|x^2) pairs - essentially
  equivalent to our gateway, but without tracking A or its growth
- arXiv:2508.07383: 4 polynomial families conjectured to cover all n equiv 1 mod 4
- arXiv:2511.07465: claims constructive proof for p equiv 1 mod 4 via affine lattice
  (preprint, convergence guarantee unclear)

### What's new (our contributions, in order of significance)
1. **Bounded-A conjecture**: max A grows by ~52 across 5 decades (199 -> 251
   from 10^6 to 10^10). No prior work identifies or discusses this phenomenon.
   If max A is bounded by a constant C, the conjecture follows from finite
   verification. Even the weaker claim max A = O(log log p) would be significant.

2. **Complete A-value distribution data**: 60.6% at A=7, 89.1% by A=11.
   This structural information about HOW primes are solved is new. Elsholtz-Tao
   count total solutions f(p), not distribution by parameter.

3. **Unified classification hierarchy**: Theorems 1-7 package cleanly, exactly
   recovering the Mordell subset. The gateway (Prop 3) provides a single
   algebraic framework for the hard cases.

4. **Constructive verification to 10^10**: Lower than 10^18, but our verification
   produces EXPLICIT decompositions (x, y, z) classified by A-value. Salez-type
   approaches confirm existence but don't classify solution structure.

5. **d | N^2 vs d | N**: The strictly weaker integrality condition is essential
   for completeness and appears not to have been explicitly noted.

### What this is NOT
- Not a proof of the full conjecture
- Not a new computational verification record (10^18 exists)
- Not claiming the gateway formulation is entirely new (2403.16047 is equivalent)

### What must be acknowledged
- arXiv:2403.16047 independently formulates the same (x, d|x^2) approach
- Our gateway parameterization by A is equivalent to their x in [ceil(p/4), ceil(p/2)]
- The novelty is tracking A and discovering its growth is bounded, not the
  algebraic identity itself

---

## Paper Structure (Revised)

### S1. Introduction (2 pages)
- State conjecture, brief history
- **Correct computational record: 10^18** (Mihnea-Dumitru 2025, using Salez algorithm)
- Position our work: "We provide a structural analysis of solutions, not a new
  verification record. Our main observation is..."
- Main result: bounded-A conjecture + classification data
- Acknowledge equivalence with 2403.16047 (x, d) formulation

### S2. Preliminaries (1 page)
- Reduction to primes (standard)
- Notation: For odd A equiv 3 (mod 4) with 4 | (p+A), set N = (p+A)/4
- Type I / Type II solution distinction (Elsholtz-Tao)

### S3. The Gateway Decomposition (2.5 pages) - CORE
- Theorem 1 (Gateway Decomposition) - the single parametric identity with d | N^2
- Remark on equivalence to 2403.16047 (x, d) formulation
- The key technical point: d | N^2 vs d | N

### S4. Existence Results (2 pages)
- Propositions for specific prime classes (~97.6% algebraic coverage)
- Remark: Propositions 2-4 together recover exactly the Mordell subset

### S5. Computational Results (2 pages)
- Verification to 10^10 with 25 A-values
- Scale progression table
- A-value distribution table
- Comparison table: our approach vs Salez modular filtering

### S6. Discussion (2.5 pages) - EXPANDED
1. **Bounded-A conjecture** (formal statement)
   - Conjecture: there exists C such that every prime admits A <= C
   - Weaker form: max A(N) = O(log log p)
   - Supporting data across 5 decades

2. **Heuristic from Elsholtz-Tao density**
   - f(p) ~ (log p)^3 explains why small A works
   - For fixed A, expected valid divisors ~ tau(N^2)/A -> infinity
   - Probability of simultaneous failure across all A decreases rapidly

3. **QR obstruction and why our approach differs**
   - Tao: QR prevents covering systems from handling Mordell subset
   - Our gateway circumvents this - uses divisor enumeration, not covering systems
   - This is why gateway can (empirically) handle cases that covering systems cannot

4. **Connection to prior work**
   - Equivalence with 2403.16047 (x, d|x^2) pairs
   - Comparison with Salez modular filtering
   - Polynomial families of 2508.07383
   - Our approach gives explicit classified decompositions

5. **The precise open problem**
   - Does N^2 always have a divisor in the target residue class for bounded A?
   - Equidistribution of divisors for structured integers (Hooley, Tenenbaum)

### S7. Conclusion (0.5 pages)

---

## Theorem Count: 5 (down from 14+)

| # | Name | Content |
|---|------|---------|
| 1 | Gateway Decomposition | The single parametric identity with d | N^2 |
| 2 | Classical (p equiv 3 mod 4) | Existence at A = 1 |
| 3 | Unified mod-3 gateway | Existence at A = 3 when B has factor not-equiv 1 mod 3 |
| 4 | NQR-7 closed forms | Existence at A = 7 for p NQR mod 7 |
| 5 | Computational coverage | 100% to 10^10 with max A = 251 |

Plus one formal **Conjecture** (bounded-A).

---

## Figures (4)

1. **Coverage hierarchy diagram**: The 4-layer classification as nested sets
2. **A-value distribution**: Bar chart for Case B QR7 primes at 10^10
3. **Max A vs limit**: Log-log showing slow growth (199 -> 251 over 5 decades)
4. **Comparison with Salez**: Our structural approach vs modular filtering

---

## Target Venue

**Experimental Mathematics** - ideal for computational-theoretical interplay.
Alternative: Mathematics of Computation.

The bounded-A conjecture as headline aligns well with Experimental Mathematics'
preference for computationally-motivated conjectures.

## Estimated Length

~12-14 pages (slightly longer than before due to expanded discussion).

---

## Pre-Writing Checklist

- [x] Literature search: Is Proposition 3 (unified mod-3 gateway) known?
      -> The A=3 special case is classical. Our unified presentation of
         the three sub-cases under one proposition may be new in packaging.
- [x] Literature search: Is Proposition 4 (NQR-7 closed forms) known?
      -> The mod-7 case is handled implicitly by Mordell's work. Our explicit
         closed forms for each NQR sub-case may be new.
- [x] Literature search: Is the gateway formulation known?
      -> YES: 2403.16047 has equivalent (x, d|x^2) formulation. Must cite.
- [x] Literature search: Is the bounded-A phenomenon noted anywhere?
      -> NO: confirmed novel. No paper mentions A-value growth or bounded parameters.
- [x] Literature search: What is the current verification record?
      -> 10^18 (Mihnea-Dumitru 2025, arXiv:2509.00128). Must cite, not compete.
- [ ] Verify the tautology proof: 4xyz = p(yz+xz+xy) from construction
- [ ] Write the d | N^2 proof carefully (the key novel claim)
- [ ] Align notation with Elsholtz & Tao (Type I/II)
- [ ] Create clean reproducibility script (single parameterized Python file)
- [ ] Generate figures from 10^10 data

---

## Bibliography (must include)

| Ref | Paper | Why |
|-----|-------|-----|
| Elsholtz-Tao 2013 | Counting solutions | f(p) bounds, Type I/II, density heuristic |
| Salez 2014 | Modular equations to 10^17 | Modular filtering comparison |
| Mihnea-Dumitru 2025 | 2509.00128 - to 10^18 | Current verification record |
| 2403.16047 | Elemental Patterns | Equivalent (x, d) formulation |
| 2508.07383 | Polynomial Families | Related polynomial approach |
| Tao 2011 | Blog: QR obstruction | Why covering systems fail |
| Mordell 1969 | Diophantine Equations | Mordell subset, reduction to primes |
| Schinzel 1956, Sierpinski 1956, Webb 1970 | Classical results | History |

---

## Writing Order

1. S3 (Gateway Decomposition) - the core, write first
2. S4 (Existence Results) - the algebraic payoff
3. S6 (Discussion) - context, bounded-A conjecture, comparisons
4. S5 (Computational Results) - the evidence
5. S1 (Introduction) - write last
6. S2, S7 - fill in
