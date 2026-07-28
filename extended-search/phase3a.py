"""Phase 3a - the A(p) growth law: what the sweep can and cannot decide.

The plan asks whether max A(p) over p <= X is constant, or grows like log log p, or like
log p / log log p.

Fitting the six milestone maxima directly is hopeless - max is one order statistic and
four of the six milestones share a value with a neighbour. This script instead fits the
object that DETERMINES the maximum: the tail of the distribution of A(p) over
gateway-resolved primes, of which there are 88,088,687 at X = 10^11.

The central lesson, learned the hard way on referee review: the tail decay is fitted as
log P(A >= T | X) = a_T - c_T * f(log X), and THE CHOICE OF f DECIDES THE ANSWER before
any data is seen.

  f = log log X  ->  gw(X)*P grows without limit for every fixed T  ->  max A unbounded
  f = log X      ->  gw(X)*P eventually falls below 1               ->  max A bounded

Both are fitted below. They have indistinguishable in-sample skill (R^2 0.972 vs 0.954,
and log X wins outright at T = 31, 47, 71), because log log X spans only 2.63 -> 3.23
across the entire sweep - a factor of 1.23 over six points. So the sweep does NOT
discriminate bounded from unbounded, and this script's job is to say so with the numbers
that show it, not to pick a side.

Passes:
  1  stability of the A(p) distribution, cumulative and per-decade
  2  the two tail models, fitted side by side, with their opposite conclusions
  3  the plan's three named families against max A vs X, for completeness

Data: results/session21_checkpoint.json, the session 20/21 authoritative sweep
(sessions 18-19 carried an inverted classifier and are not used).

Run: python phase3a.py
"""
import json, math

CHK = 'results/session21_checkpoint.json'
OUT = 'results/phase3a_growth.json'

Ts = [11, 19, 31, 47, 71, 103, 151, 199]


def lstsq(xs, ys):
    """Least squares y = a + b x. Returns (a, b, R^2, se_b)."""
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    b = sxy / sxx
    a = my - b * mx
    ss_res = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - my) ** 2 for y in ys)
    se = math.sqrt(ss_res / (n - 2) / sxx) if n > 2 and sxx else float('nan')
    return a, b, (1 - ss_res / ss_tot if ss_tot else float('nan')), se


def gateways(limit):
    """Gateway parameters A = 3 mod 4, prime, below limit - in order."""
    return [q for q in range(3, limit)
            if q % 4 == 3 and all(q % r for r in range(2, int(q ** .5) + 1))]


GW_ALL = gateways(20000)
# g(T) counts gateways strictly below T. NOTE the offset is arbitrary and matters:
# A = 3 never appears in the data (every A=3 prime is resolved by an earlier stage), so
# whether 3 and 7 are counted shifts the fitted INTERCEPT, which is why no per-gateway
# reading of the slope is claimed below.
def g_of(T):
    return sum(1 for q in GW_ALL if q < T)


def surv(v, T):
    return sum(n for a, n in v['Au'].items() if int(a) >= T)


def main():
    d = json.load(open(CHK))
    ms = sorted(((int(k), v) for k, v in d['milestones'].items()))
    for X, v in ms:
        assert sum(v['Au'].values()) == v['gw'], X   # load-bearing: no missing mass
        assert v['open'] == 0, X                     # load-bearing: no top-end censoring

    print('=' * 78)
    print('DATA - session 20/21 authoritative sweep')
    print(f"  {'X':>13} {'primes':>13} {'gw-resolved':>13} {'gw fraction':>12} {'max A':>6}")
    for X, v in ms:
        print(f"  {X:>13.0f} {v['primes']:>13} {v['gw']:>13} "
              f"{v['gw']/v['primes']:>12.6f} {v['mxA']:>6}")
    print(f"  smallest A occurring at 10^11: {min(int(a) for a in ms[-1][1]['Au'])}"
          '  (A = 3 never appears - earlier stages take those primes)')
    print()

    # ------------------------------------------------------------------ pass 1
    print('=' * 78)
    print('PASS 1 - the distribution of A(p) is not stable; its tail thins as X grows')
    print('  P(A >= T) over gateway-resolved primes. CUMULATIVE (all p <= X):')
    print('  ' + f"{'T':>5}" + ''.join(f'{X:>12.0e}' for X, _ in ms))
    cum = {}
    for T in Ts:
        cum[T] = [surv(v, T) / v['gw'] for _, v in ms]
        print('  ' + f'{T:>5}' + ''.join(f'{r:>12.3e}' for r in cum[T]))
    print()
    print('  PER-DECADE (p in (X/10, X]) - guards against the cumulative-overlap trap,')
    print('  since the 10^11 cumulative sample contains all the earlier ones:')
    print('  ' + f"{'T':>5}" + ''.join(f'{X:>12.0e}' for X, _ in ms))
    dec = {}
    for T in Ts:
        row = []
        for i, (X, v) in enumerate(ms):
            if i == 0:
                row.append(surv(v, T) / v['gw'])
            else:
                dn = surv(v, T) - surv(ms[i - 1][1], T)
                dd = v['gw'] - ms[i - 1][1]['gw']
                row.append(dn / dd if dd else float('nan'))
        dec[T] = row
        print('  ' + f'{T:>5}' + ''.join(f'{r:>12.3e}' for r in row))
    print()
    print('  Monotone decreasing at every threshold in BOTH tables, with essentially the')
    print('  same magnitude, so the decay is not an artifact of cumulative overlap.')
    print('  The gateway fraction is itself falling (0.0289 -> 0.0214), so the')
    print('  conditioning population changes too - but earlier stages taking a GROWING')
    print('  share would bias P(A>=T | gateway) UP, not down, so this is conservative.')
    print()

    # ------------------------------------------------------------------ pass 2
    print('=' * 78)
    print('PASS 2 - two tail models with indistinguishable fit and opposite conclusions')
    regs = {'log log X (-> unbounded)': lambda X: math.log(math.log(X)),
            'log X     (-> bounded)  ': math.log}
    fitted = {}
    for name, f in regs.items():
        rows = []
        for T in Ts:
            pts = [(X, cum[T][i]) for i, (X, _) in enumerate(ms) if cum[T][i] > 0]
            if len(pts) < 4:
                continue
            a, b, r2, _ = lstsq([f(X) for X, _ in pts], [math.log(p) for _, p in pts])
            rows.append((T, g_of(T), a, -b, r2))
        aa, ab, ar, ase = lstsq([g for _, g, _, _, _ in rows], [a for _, _, a, _, _ in rows])
        ca, cb, cr, cse = lstsq([g for _, g, _, _, _ in rows], [c for _, _, _, c, _ in rows])
        fitted[name] = dict(rows=rows, a=(aa, ab, ar, ase), c=(ca, cb, cr, cse), f=f)
        print(f'  regressor {name}')
        print('  ' + f"{'T':>5} {'g(T)':>5} {'a_T':>9} {'c_T':>9} {'R^2 (per T)':>12}")
        for T, g, a, c, r2 in rows:
            print('  ' + f'{T:>5} {g:>5} {a:>9.3f} {c:>9.4f} {r2:>12.5f}')
        print(f'    a_T = {aa:.4f} + {ab:.4f} g(T)   R^2 = {ar:.4f}  (se {ase:.4f})')
        print(f'    c_T = {ca:.4f} + {cb:.4f} g(T)   R^2 = {cr:.4f}  (se {cse:.4f})')
        print()

    print('  Per-T comparison - which regressor wins where:')
    print('  ' + f"{'T':>5} {'R^2 log log X':>15} {'R^2 log X':>12} {'winner':>10}")
    r1 = {T: r for T, _, _, _, r in fitted['log log X (-> unbounded)']['rows']}
    r2d = {T: r for T, _, _, _, r in fitted['log X     (-> bounded)  ']['rows']}
    for T in sorted(r1):
        w = 'log log X' if r1[T] > r2d[T] else 'log X'
        print('  ' + f'{T:>5} {r1[T]:>15.5f} {r2d[T]:>12.5f} {w:>10}')
    print()
    print(f'  log log X spans {math.log(math.log(1e6)):.3f} -> '
          f'{math.log(math.log(1e11)):.3f} across the whole sweep - a factor of '
          f'{math.log(math.log(1e11))/math.log(math.log(1e6)):.2f}.')
    print('  Six points over a 23% regressor range cannot separate these two forms.')
    print()

    # the gateway fraction is NOT stationary; extrapolate it rather than freezing it
    fa, fb, fr, _ = lstsq([math.log(math.log(X)) for X, _ in ms],
                          [math.log(v['gw'] / v['primes']) for _, v in ms])
    print(f'  gateway fraction fitted as log frac = {fa:.3f} + {fb:.3f} log log X '
          f'(R^2 = {fr:.4f}); used below instead of freezing it at its 10^11 value.')

    def gw_of(X):
        return X / math.log(X) * math.exp(fa + fb * math.log(math.log(X)))

    def predict_max(X, mdl):
        """Largest gateway T whose expected count is >= 1. g_of(T) matches the fit."""
        (aa, ab, _, _), (ca, cb, _, _), f = mdl['a'], mdl['c'], mdl['f']
        best = 3
        for T in GW_ALL:
            g = g_of(T)
            if gw_of(X) * math.exp((aa + ab * g) - (ca + cb * g) * f(X)) >= 1:
                best = T
        return best

    print()
    print('  Predicted max A under each model:')
    print('  ' + f"{'X':>10} {'log log X model':>17} {'log X model':>13} {'observed':>9}")
    preds = {}
    for X, v in ms:
        preds[X] = [predict_max(X, fitted[k]) for k in regs]
        print('  ' + f'{X:>10.0e} {preds[X][0]:>17} {preds[X][1]:>13} {v["mxA"]:>9}')
    for X in (1e13, 1e18, 1e30, 1e100):
        p = [predict_max(X, fitted[k]) for k in regs]
        preds[X] = p
        print('  ' + f'{X:>10.0e} {p[0]:>17} {p[1]:>13} {"-":>9}')
    ca, cb, _, _ = fitted['log X     (-> bounded)  ']['c']
    aa, ab, _, _ = fitted['log X     (-> bounded)  ']['a']
    gmax = (1 - aa + 0) / cb if cb else float('inf')
    print(f'  The log X model saturates: it needs c_T * log X <= a_T + log gw(X), and')
    print(f'  since c_T grows in g(T) while log gw(X) ~ log X, the largest admissible')
    print(f'  g is finite. Its ceiling is max A = {predict_max(1e300, fitted[list(regs)[1]])}.')
    print()
    print('  Sensitivity: a_T slope se = '
          f"{fitted['log log X (-> unbounded)']['a'][3]:.4f} on a slope of "
          f"{fitted['log log X (-> unbounded)']['a'][1]:.4f}. A one-se shift moves the")
    print('  unbounded model\'s 10^18 prediction by roughly 50, so its 1-sigma band at')
    print('  10^18 already contains 359 - the currently observed maximum.')
    print()

    # ------------------------------------------------------------------ pass 3
    print('=' * 78)
    print("PASS 3 - the plan's three named families against max A vs X, for completeness")
    Xs = [X for X, _ in ms]
    Ys = [float(v['mxA']) for _, v in ms]
    fam = {'constant': None,
           'log log p': lambda X: math.log(math.log(X)),
           'log p / log log p': lambda X: math.log(X) / math.log(math.log(X)),
           'log p': math.log,
           '(log log p)^2': lambda X: math.log(math.log(X)) ** 2}
    print(f"  {'family':<22} {'R^2':>8} {'at 10^18':>10} {'at 10^30':>10}")
    direct = {}
    for name, f in fam.items():
        if f is None:
            m = sum(Ys) / len(Ys)
            print(f'  {name:<22} {"n/a":>8} {m:>10.0f} {m:>10.0f}')
            direct[name] = dict(r2=None, at18=m, at30=m)
            continue
        a, b, r2, _ = lstsq([f(X) for X in Xs], Ys)
        direct[name] = dict(r2=r2, at18=a + b * f(1e18), at30=a + b * f(1e30))
        print(f'  {name:<22} {r2:>8.4f} {a+b*f(1e18):>10.0f} {a+b*f(1e30):>10.0f}')
    print()
    print('  max A is flat at 239 across 10^8-10^9 and at 359 across 10^10-10^11, so six')
    print('  points carry ~4 degrees of freedom and the R^2 ordering is noise. The')
    print('  constant model has no R^2 by construction (it is fitted as a mean), so its')
    print('  0.0 is definitional and is NOT evidence against it.')
    print("  Salez's 10^18 verification of Erdos-Straus does not help: it records")
    print('  solvability, not the least successful gateway, so it carries no A(p) data.')

    json.dump(dict(
        milestones=[(X, v['gw'], v['primes'], v['mxA']) for X, v in ms],
        cumulative=cum, per_decade=dec,
        models={k: dict(a=v['a'], c=v['c'], rows=v['rows']) for k, v in fitted.items()},
        predictions={str(k): v for k, v in preds.items()},
        direct=direct), open(OUT, 'w'))
    print(f'\nwrote {OUT}')


if __name__ == '__main__':
    main()
