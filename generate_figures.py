#!/usr/bin/env python3
"""Generate all paper figures from collected statistics."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'figure.dpi': 200,
    'savefig.dpi': 200,
    'savefig.bbox': 'tight',
})

with open("stats_1M.json") as f:
    stats = json.load(f)

# ============================================================================
# FIGURE 1: Coverage Hierarchy — Nested Pie / Sunburst
# ============================================================================

def fig1_coverage():
    cov = stats["coverage"]
    total = stats["total"]

    labels = [
        r'$p \equiv 3\;(\mathrm{mod}\;4)$',
        r'$p \equiv 5\;(\mathrm{mod}\;8)$',
        r'$p \equiv 17\;(\mathrm{mod}\;24)$',
        'Case A',
        'Case B\nNQR mod 7',
        'Case B\nQR mod 7',
    ]
    values = [
        cov["3mod4"], cov["5mod8"], cov["17mod24"],
        cov["caseA"], cov["nqr7"], cov["qr7"]
    ]
    colors = [
        '#2196F3',  # blue
        '#4CAF50',  # green
        '#FF9800',  # orange
        '#9C27B0',  # purple
        '#F44336',  # red
        '#E91E63',  # pink
    ]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 7),
                                     gridspec_kw={'wspace': 0.3})

    # Left: Pie chart — start at 220° so small slices end up on the lower-right
    wedges, texts, autotexts = ax1.pie(
        values, labels=None, colors=colors, autopct='',
        startangle=220, counterclock=False,
        wedgeprops={'linewidth': 1.5, 'edgecolor': 'white'},
        radius=1.0
    )
    # Inside labels for big slices only
    for i, (w, v) in enumerate(zip(wedges, values)):
        pct = 100 * v / total
        mid_angle = (w.theta2 + w.theta1) / 2
        if pct > 8:
            r = 0.6
            x = r * np.cos(np.radians(mid_angle))
            y = r * np.sin(np.radians(mid_angle))
            ax1.text(x, y, f'{pct:.1f}%', ha='center', va='center',
                    fontsize=11, color='white', fontweight='bold')

    ax1.set_title('Prime Classification\n(all primes to $10^6$)', pad=12, fontsize=13)

    # Legend below the left pie
    legend_labels = [f'{labels[i]}  ({100*values[i]/total:.1f}%)'
                     for i in range(len(labels))]
    ax1.legend(wedges, legend_labels, loc='upper center',
              bbox_to_anchor=(0.5, -0.05), fontsize=8.5, frameon=False,
              ncol=2, columnspacing=1.0)

    # Right: Zoom into the ~5% "hard" primes
    hard_labels_short = ['Case A', 'NQR mod 7', 'QR mod 7']
    hard_refs = ['(Prop. 4.2c)', '(Prop. 4.3)', '(Thm. 5.1)']
    hard_values = [cov["caseA"], cov["nqr7"], cov["qr7"]]
    hard_colors = ['#9C27B0', '#F44336', '#E91E63']

    wedges2, _, _ = ax2.pie(
        hard_values, labels=None, colors=hard_colors, autopct='',
        startangle=90, counterclock=False,
        wedgeprops={'linewidth': 1.5, 'edgecolor': 'white'},
        radius=1.0
    )
    for i, (w, v) in enumerate(zip(wedges2, hard_values)):
        pct = 100 * v / sum(hard_values)
        angle = (w.theta2 + w.theta1) / 2
        x = 0.55 * np.cos(np.radians(angle))
        y = 0.55 * np.sin(np.radians(angle))
        ax2.text(x, y, f'{pct:.0f}%\n({v:,})', ha='center', va='center',
                fontsize=9, color='white', fontweight='bold')

    ax2.set_title('Zoom: $p \\equiv 1$ (mod 24)\nCase B primes',
                  pad=12, fontsize=13)

    legend2 = [f'{hard_labels_short[i]}  {hard_refs[i]}  ({hard_values[i]:,})'
               for i in range(len(hard_labels_short))]
    ax2.legend(wedges2, legend2, loc='upper center',
              bbox_to_anchor=(0.5, -0.05), fontsize=8.5, frameon=False)

    plt.savefig('fig1_coverage.png', bbox_inches='tight')
    plt.savefig('fig1_coverage.pdf', bbox_inches='tight')
    print("Figure 1 saved: fig1_coverage.png/pdf")
    plt.close()


# ============================================================================
# FIGURE 2: A-value Distribution (bar chart)
# ============================================================================

def fig2_a_distribution():
    A_dist = {int(k): v for k, v in stats["A_dist"].items()}
    total_qr7 = stats["coverage"]["qr7"]

    # Sort by A value
    A_vals = sorted(A_dist.keys())
    counts = [A_dist[a] for a in A_vals]
    pcts = [100 * c / total_qr7 for c in counts]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), height_ratios=[3, 1])

    # Top: bar chart
    bars = ax1.bar(range(len(A_vals)), pcts, color='#2196F3', edgecolor='#1565C0',
                   linewidth=0.5)
    ax1.set_xticks(range(len(A_vals)))
    ax1.set_xticklabels([str(a) for a in A_vals], rotation=45, ha='right', fontsize=8)
    ax1.set_ylabel('% of Case B QR7 primes solved')
    ax1.set_title(f'Distribution of first-working $A$ value\n'
                  f'({total_qr7:,} Case B QR7 primes to $10^6$)')
    ax1.set_xlim(-0.5, len(A_vals) - 0.5)

    # Add value labels on top of bars > 1%
    for i, (bar, pct, cnt) in enumerate(zip(bars, pcts, counts)):
        if pct > 1:
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    f'{pct:.1f}%', ha='center', va='bottom', fontsize=8)

    # Bottom: cumulative
    cum = np.cumsum(pcts)
    ax2.fill_between(range(len(A_vals)), cum, alpha=0.3, color='#2196F3')
    ax2.plot(range(len(A_vals)), cum, 'o-', color='#1565C0', markersize=4, linewidth=1.5)
    ax2.set_xticks(range(len(A_vals)))
    ax2.set_xticklabels([str(a) for a in A_vals], rotation=45, ha='right', fontsize=8)
    ax2.set_ylabel('Cumulative %')
    ax2.set_ylim(0, 105)
    ax2.axhline(y=100, color='gray', linestyle='--', linewidth=0.5)
    ax2.set_xlim(-0.5, len(A_vals) - 0.5)
    ax2.set_xlabel('$A$ value (prime, $A \\equiv 3\\;(\\mathrm{mod}\\;4)$)')

    # Annotate key points
    for i, a in enumerate(A_vals):
        if cum[i] >= 99.9 and (i == 0 or cum[i-1] < 99.9):
            ax2.annotate(f'99.9% at $A={a}$', xy=(i, cum[i]),
                        xytext=(i+2, 85), fontsize=9,
                        arrowprops=dict(arrowstyle='->', color='red'),
                        color='red')
            break

    plt.tight_layout()
    plt.savefig('fig2_A_distribution.png')
    plt.savefig('fig2_A_distribution.pdf')
    print("Figure 2 saved: fig2_A_distribution.png/pdf")
    plt.close()


# ============================================================================
# FIGURE 3: Max A Stabilisation (across scales)
# ============================================================================

def fig3_max_A_stabilisation():
    # Data from sessions 14-17 (10^3..10^9) and the corrected session-20
    # re-verification (10^10); see REFEREE_FINDINGS_2026-07-05.md
    scale_data = {
        1_000: 7,
        2_000: 7,
        5_000: 23,
        10_000: 31,
        20_000: 47,
        50_000: 47,
        100_000: 59,
        200_000: 59,
        500_000: 79,
        1_000_000: 79,
        10_000_000: 167,
        100_000_000: 239,
        1_000_000_000: 239,
        10_000_000_000: 359,
        100_000_000_000: 359,
    }
    # Override with our progressive data where available
    prog = stats.get("max_A_progressive", {})
    for k, v in prog.items():
        scale_data[int(k)] = v

    limits = sorted(scale_data.keys())
    max_As = [scale_data[l] for l in limits]

    fig, ax = plt.subplots(figsize=(9, 5))

    ax.semilogx(limits, max_As, 'o-', color='#E91E63', markersize=7,
                linewidth=2, markeredgecolor='#880E4F', markeredgewidth=1)

    # Reference lines at the plateau and the current maximum
    ax.axhline(y=239, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.text(2e9, 239, '$A = 239$', va='bottom', ha='right', fontsize=9,
            color='gray', style='italic')
    ax.axhline(y=359, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.text(2e11, 359, '$A = 359$', va='bottom', ha='right', fontsize=9,
            color='gray', style='italic')

    # Logarithmic growth guide, then the flat 10^10 -> 10^11 step
    ax.text(3e7, 30, 'Roughly linear in $\\log N$', fontsize=9,
            color='#880E4F', alpha=0.8)
    ax.annotate('no growth\n$10^{10}\\to10^{11}$', xy=(1e11, 359),
                xytext=(4e9, 300), fontsize=8, color='#880E4F',
                ha='center', arrowprops=dict(arrowstyle='->', color='#880E4F',
                                             alpha=0.7))

    ax.set_xlabel('Verification limit')
    ax.set_ylabel('Maximum $A$ needed')
    ax.set_title('Growth of max $A$ across scales')
    ax.set_ylim(0, 400)
    ax.set_xlim(500, 4e11)

    # Custom x-ticks
    xticks = [1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11]
    ax.set_xticks(xticks)
    ax.set_xticklabels(['$10^3$', '$10^4$', '$10^5$', '$10^6$',
                        '$10^7$', '$10^8$', '$10^9$', '$10^{10}$', '$10^{11}$'])
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('fig3_max_A.png')
    plt.savefig('fig3_max_A.pdf')
    print("Figure 3 saved: fig3_max_A.png/pdf")
    plt.close()


# ============================================================================
# FIGURE 4: d/N ratio — showing the N^2 extension
# ============================================================================

def fig4_d_over_N():
    data = stats["d_over_N"]
    ratios = [x["ratio"] for x in data]
    primes = [x["p"] for x in data]
    divides_N = [x["d_divides_N"] for x in data]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left: scatter of d/N ratio vs p
    colors = ['#2196F3' if dn else '#E91E63' for dn in divides_N]
    sizes = [3 if dn else 12 for dn in divides_N]

    # Plot d|N points first (smaller, blue), then d∤N points on top (larger, pink)
    for is_div in [True, False]:
        idx = [i for i in range(len(data)) if divides_N[i] == is_div]
        ps = [primes[i] for i in idx]
        rs = [ratios[i] for i in idx]
        c = '#2196F3' if is_div else '#E91E63'
        s = 3 if is_div else 15
        label = f'$d \\mid N$ ({len(idx):,})' if is_div else f'$d \\nmid N$ ({len(idx):,})'
        alpha = 0.3 if is_div else 0.8
        ax1.scatter(ps, rs, c=c, s=s, alpha=alpha, label=label, zorder=3 if not is_div else 2)

    ax1.axhline(y=1.0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax1.text(50000, 1.05, '$d/N = 1$ (boundary)', fontsize=9, color='gray')
    ax1.set_xlabel('Prime $p$')
    ax1.set_ylabel('$d / N$')
    ax1.set_title('Gateway divisor ratio $d/N$\n(Case B QR7 primes to $10^6$)')
    ax1.legend(fontsize=9, loc='upper right')
    ax1.set_yscale('log')
    ax1.set_ylim(0.001, 100)

    # Right: histogram of log(d/N)
    log_ratios = [np.log10(r) for r in ratios if r > 0]
    log_ratios_divN = [np.log10(ratios[i]) for i in range(len(data))
                       if divides_N[i] and ratios[i] > 0]
    log_ratios_notN = [np.log10(ratios[i]) for i in range(len(data))
                       if not divides_N[i] and ratios[i] > 0]

    bins = np.linspace(-3, 2, 40)
    ax2.hist(log_ratios_divN, bins=bins, color='#2196F3', alpha=0.6,
             label=f'$d \\mid N$ ({len(log_ratios_divN):,})', edgecolor='white')
    ax2.hist(log_ratios_notN, bins=bins, color='#E91E63', alpha=0.8,
             label=f'$d \\nmid N$ ({len(log_ratios_notN):,})', edgecolor='white')

    ax2.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax2.text(0.05, ax2.get_ylim()[1]*0.9 if ax2.get_ylim()[1] > 0 else 100,
             '$d = N$', fontsize=9, color='gray')
    ax2.set_xlabel('$\\log_{10}(d/N)$')
    ax2.set_ylabel('Count')
    ax2.set_title('Distribution of $d/N$ ratios')
    ax2.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig('fig4_d_over_N.png')
    plt.savefig('fig4_d_over_N.pdf')
    print("Figure 4 saved: fig4_d_over_N.png/pdf")
    plt.close()


# ============================================================================
# FIGURE 5 (BONUS): The Gateway Formula — visual diagram
# ============================================================================

def fig5_gateway_diagram():
    """Conceptual diagram of the gateway decomposition."""
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    ax.axis('off')

    # Title
    ax.text(5, 4.5, 'The Gateway Decomposition', ha='center', fontsize=16,
            fontweight='bold', style='italic')

    # Main equation
    ax.text(5, 3.8,
            r'$\frac{4}{p} = \frac{1}{x} + \frac{1}{y} + \frac{1}{z}$',
            ha='center', fontsize=18, math_fontfamily='cm')

    # Arrow down
    ax.annotate('', xy=(5, 3.0), xytext=(5, 3.5),
                arrowprops=dict(arrowstyle='->', color='#333', lw=2))

    # Three boxes for x, y, z
    box_style = dict(boxstyle='round,pad=0.4', facecolor='#E3F2FD',
                     edgecolor='#1565C0', linewidth=1.5)
    ax.text(2, 2.5, r'$x = \frac{p + A}{4}$', ha='center', fontsize=14,
            bbox=box_style, math_fontfamily='cm')
    ax.text(5, 2.5, r'$y = \frac{pN + d}{A}$', ha='center', fontsize=14,
            bbox=box_style, math_fontfamily='cm')
    ax.text(8, 2.5, r'$z = \frac{pNy}{d}$', ha='center', fontsize=14,
            bbox=box_style, math_fontfamily='cm')

    # Conditions
    cond_style = dict(boxstyle='round,pad=0.3', facecolor='#FCE4EC',
                      edgecolor='#C62828', linewidth=1.2)
    ax.text(2, 1.2, r'$4 \mid (p + A)$', ha='center', fontsize=11,
            bbox=cond_style, math_fontfamily='cm')
    ax.text(5, 1.2, r'$A \mid (pN + d)$', ha='center', fontsize=11,
            bbox=cond_style, math_fontfamily='cm')
    ax.text(8, 1.2, r'$d \mid N^2$', ha='center', fontsize=11,
            bbox=cond_style, math_fontfamily='cm')

    # Labels
    ax.text(2, 0.6, 'Divisibility', ha='center', fontsize=9, color='#C62828')
    ax.text(5, 0.6, 'Residue', ha='center', fontsize=9, color='#C62828')
    ax.text(8, 0.6, 'Gateway', ha='center', fontsize=9, color='#C62828')

    # Connecting arrows
    for cx in [2, 5, 8]:
        ax.annotate('', xy=(cx, 1.55), xytext=(cx, 2.1),
                    arrowprops=dict(arrowstyle='->', color='#999', lw=1))

    # Where N = ...
    ax.text(5, 0.15, r'where $N = (p+A)/4$,  $A$ prime,  $A \equiv 3\;(\mathrm{mod}\;4)$',
            ha='center', fontsize=10, color='#555')

    plt.tight_layout()
    plt.savefig('fig5_gateway_diagram.png')
    plt.savefig('fig5_gateway_diagram.pdf')
    print("Figure 5 saved: fig5_gateway_diagram.png/pdf")
    plt.close()


# Run all
fig1_coverage()
fig2_a_distribution()
fig3_max_A_stabilisation()
fig4_d_over_N()
fig5_gateway_diagram()
print("\nAll figures generated!")
