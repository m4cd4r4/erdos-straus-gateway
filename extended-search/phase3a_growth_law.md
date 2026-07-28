# Phase 3a - the A(p) growth law: what the sweep can and cannot decide

**Run 2026-07-25. Script `phase3a.py`, records `results/phase3a_growth.json`, console log
`results/phase3a_run.txt`. Data: `results/session21_checkpoint.json`, the session 20/21
authoritative sweep (sessions 18-19 carried an inverted classifier and are not used).**

**Refereed 2026-07-25 by a fresh-context adversarial reviewer: PARTIALLY REFUTED. The
first draft's headline - "A(p) is almost certainly NOT bounded, and the growth is
Stormer-shaped" - is WITHDRAWN in full. It was an artifact of the fitted functional form,
not a finding. See "What the first draft got wrong". The referee also found this project's
first genuine code defect, an off-by-one; it is fixed.**

## Verdict

**The sweep cannot discriminate bounded from unbounded A(p). Neither the paper's
bounded-A headline nor its negation is supported by this data.** What *is* established is
narrower and still useful: the distribution of `A(p)` is not stationary, its tail thins as
`X` grows, and therefore **the flatness of `max A` since `10^10` is not by itself evidence
of a ceiling** - which is how 22 sessions of sweeping had been reading it.

## Finding 1 - the distribution of A(p) is not stable; its tail thins as X grows

`P(A(p) >= T)` over gateway-resolved primes, **cumulative** (all `p <= X`):

| T | 10^6 | 10^7 | 10^8 | 10^9 | 10^10 | 10^11 |
|---|---|---|---|---|---|---|
| 11 | 4.359e-1 | 3.965e-1 | 3.708e-1 | 3.476e-1 | 3.276e-1 | 3.109e-1 |
| 31 | 3.482e-2 | 2.897e-2 | 2.257e-2 | 1.794e-2 | 1.434e-2 | 1.147e-2 |
| 71 | 3.085e-3 | 1.942e-3 | 1.210e-3 | 8.248e-4 | 5.368e-4 | 3.469e-4 |

and **per decade** (`p` in `(X/10, X]`), which is the version that matters, since the
cumulative `10^11` sample contains every earlier sample and successive points are
therefore massively correlated:

| T | 10^6 | 10^7 | 10^8 | 10^9 | 10^10 | 10^11 |
|---|---|---|---|---|---|---|
| 11 | 4.359e-1 | 3.908e-1 | 3.672e-1 | 3.444e-1 | 3.250e-1 | 3.087e-1 |
| 31 | 3.482e-2 | 2.813e-2 | 2.167e-2 | 1.730e-2 | 1.386e-2 | 1.110e-2 |
| 71 | 3.085e-3 | 1.778e-3 | 1.108e-3 | 7.719e-4 | 4.981e-4 | 3.219e-4 |

Monotone decreasing at every threshold in both tables, with the same magnitude. The decay
is real and is not a cumulative-overlap artifact.

Two confounds, both checked:

- **The denominator moves too.** The gateway-resolved fraction of primes falls from
  0.028905 to 0.021391 across the sweep. But earlier resolution stages taking a *growing*
  share of primes would bias `P(A >= T | gateway)` **up**, not down, so the conditioning
  works against the observed decay. Conservative.
- **The decay survives per-decade refitting**, so it is not an artifact of the correlation
  between cumulative milestones.

**Consequence.** Any extreme-value analysis that treats `A(p)` as having a fixed limiting
distribution will mis-price the maximum. A stationary tail fitted at `10^11` under-predicts
the observed maxima at every smaller milestone. (The first draft attributed all 27% of that
error to non-stationarity; the referee correctly notes it under-predicts by 15% at `10^11`
itself, where stationarity holds by construction, so most of the gap is extreme-value
variance and only some is non-stationarity.)

## Finding 2 - the tail decays as a power of a slowly-growing function of X, with an exponent rising in T

Fitting `log P(A >= T | X) = a_T - c_T f(log X)`, the exponent `c_T` rises monotonically
with `T` under either choice of `f`, and is close to linear in `g(T) = #{gateways < T}`:

| regressor | fitted `c_T` | R^2 |
|---|---|---|
| `log log X` | `-0.0124 + 0.3431 g(T)` | 0.9718 |
| `log X` | `0.0118 + 0.0157 g(T)` | 0.9542 |

That much is real: **the more gateways a prime must fail, the faster its probability
decays in X**, which is qualitatively what the smoothness density predicts.

**Three things the first draft claimed here are withdrawn:**

- *"a power of `log X` specifically"* - unidentified. The `log X` regressor beats
  `log log X` at `T` = 31, 47, 71 and 103, and loses at 11, 19, 151, 199. Neither form is
  selected by the data.
- *"0.343 per gateway"* - not identified as a per-gateway cost. `c_T` fits `T` at
  R^2 = 0.9668 and `T/log T` at 0.9704, against 0.9718 for `g(T)`: any monotone regressor
  fits eight monotone points. Worse, **A = 3 never occurs in the data** (the smallest `A`
  present at `10^11` is 7), so `g`'s offset is arbitrary - dropping `A = 3` moves the
  intercept from `-0.012` to `0.331`, and dropping 3 and 7 moves it to `0.674`. The
  near-zero intercept that made "per gateway" readable was an artifact of that offset.
- *"0.343 vs the 0.5 of independence, therefore positive correlation"* - numerology, given
  the above. Withdrawn.

## Finding 3 - the two models are indistinguishable in-sample and opposite out of it

This is the actual result, and it is a negative one.

The conclusion is decided by the choice of regressor *before any data is seen*:

- `f = log log X`: `gw(X) * P(A >= T | X)` grows without limit at every fixed `T`, so
  **max A is unbounded**. No dataset could return "bounded" from this parameterisation.
- `f = log X`: the count eventually falls below 1, so **max A saturates**.

Both fitted, using the same pipeline, with the gateway fraction extrapolated rather than
frozen (it fits `log frac = -2.219 - 0.502 log log X`, R^2 = 0.9974 - an exponent
strikingly close to the `(log X)^{-1/2}` of the smoothness density):

| X | `log log X` model | `log X` model | observed |
|---|---|---|---|
| 10^6 | 127 | 127 | 79 |
| 10^7 | 151 | 151 | 167 |
| 10^8 | 167 | 167 | 239 |
| 10^9 | 199 | 199 | 239 |
| 10^10 | 223 | 223 | 359 |
| 10^11 | 239 | 239 | 359 |
| 10^13 | 283 | 271 | - |
| 10^18 | 419 | 359 | - |
| 10^30 | 647 | 463 | - |
| 10^100 | 2027 | 599 | - |

**The two models agree exactly at every point where data exists**, and diverge only where
none does. The `log X` model's ceiling is `max A = 643`, forever. The `log log X` model
grows without limit. Nothing in the sweep chooses between them, because `log log X` spans
only `2.626 -> 3.232` across the entire range - a factor of 1.23 over six points.

The unbounded model's own 1-sigma band at `10^18` (`a_T` slope 0.4795, se 0.0900) already
contains 359, the currently observed maximum. Both models also under-predict at `10^10`
and `10^11` by about a third, so neither is trustworthy at the level even in-sample.

## What this means for the paper

Weaker than the first draft claimed, and still worth writing down:

1. **Bounded-A is not refuted.** A saturating model fits the tail as well as a growing one
   and holds `359` through `10^18`.
2. **Bounded-A is also not supported by flatness.** "Flat since `10^10`" is exactly what
   both models predict over one decade, including the unbounded one. The flatness
   argument, which is the empirical case in the paper, carries no weight on its own.
3. **The discriminating experiment is small and identified.** The two models separate by
   only ~12 at `10^13` and 60 at `10^18`, so extending the sweep buys less than the first
   draft implied - but a single new record above 359 before `10^13` would kill the
   saturating model outright, and the sweep code is reusable and checkpointed.

This sits beside Phase 4.2 rather than against it: 4.2 showed the constructor has no lever
cheaper than Q_A-smoothness, so large `A(p)` cannot be *forced* cheaply. Phase 3a shows the
sweep cannot tell whether it nevertheless *drifts* upward. Both are statements about what
is not known.

## What the first draft got wrong

Recorded deliberately; this project has shipped an inverted classifier and partially
retracted two prior analyses, all of which passed self-review.

1. **Headline "A(p) is almost certainly NOT bounded" and "rules out constant and
   log log p".** Circular: unboundedness is guaranteed by the `log log X` parameterisation
   for any finite fitted coefficients. The referee constructed the bounded alternative,
   with equal in-sample skill and a ceiling of 643. Fully withdrawn.
2. **"Rules out constant"** leaned on the constant model's `R^2 = 0.0` in Pass 3, which is
   definitional - it is fitted as a mean - not evidence. Pass 3 now says so.
3. **"0.343 per gateway"** and the comparison to 0.5: unidentified, and the near-zero
   intercept was an artifact of counting `A = 3`, which never occurs in the data.
4. **"The shape is robust even though the level isn't"** - having it both ways. The shape
   was assumed by the functional form; the level was the only thing the data touched, and
   it is wrong by up to 37% in-sample.
5. **Froze the gateway fraction** at its `10^11` value while making non-stationarity of the
   numerator the headline finding. It falls ~5.9% per decade; carried to `10^18` that
   over-stated `gw` by ~1.5x and biased every extrapolation upward. Now fitted.
6. **CODE DEFECT (the project's first).** `predict_max` incremented `g` before use, so it
   evaluated the model at `g(T)+1` while the fit had been given `g(T)`. Every prediction
   was low by one gateway step. Fixed; `g_of(T)` is now shared by both.
7. **Selective table** - displayed 7 of the 8 fitted `T` values, omitting `T = 103`, the
   one point that breaks the monotone rise of `c_T`. All 8 are now shown.

## Limits

- Six milestones, `10^6` to `10^11`, one sweep. Every figure past `10^11` is a fitted
  surface talking.
- The population is gateway-resolved primes, already conditioned on failing thm1, thm2,
  thm4 and caseA. That conditioning is not disentangled anywhere in this file.
- Expected-count-equals-one is a crude estimator for a maximum; it ignores the variance of
  the extreme, which is large exactly where the count is O(1).
- Nothing here bears on the Section 5 character-sum problem, which remains open.

## Next step

Extend the sweep to `10^12-10^13`. It is the only cheap experiment that moves this, and
the outcome is sharply defined: **any prime with `A(p) > 359` below `10^13` kills the
saturating model**, while another two decades of flatness at 359 is mild evidence for it.
Expect the models to differ by only ~12 at `10^13`, so this is a record hunt, not a
curve-fitting exercise.
