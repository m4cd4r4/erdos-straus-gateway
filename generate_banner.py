#!/usr/bin/env python3
"""Generate a composite banner image for the README."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

with open("stats_1M.json") as f:
    stats = json.load(f)

fig = plt.figure(figsize=(16, 5.5))
fig.patch.set_facecolor('#0d1117')

gs = gridspec.GridSpec(1, 3, width_ratios=[1.1, 1.2, 1.0], wspace=0.32)

# Shared style
GOLD = '#f0c040'
CYAN = '#58a6ff'
PINK = '#f778ba'
WHITE = '#e6edf3'
DIMMED = '#8b949e'
BG = '#0d1117'
CARD = '#161b22'

# ============================================================================
# PANEL 1: Coverage pie - simplified single pie
# ============================================================================
ax1 = fig.add_subplot(gs[0])
ax1.set_facecolor(BG)

cov = stats["coverage"]
total = stats["total"]

labels = [
    r'$p \equiv 3\;(4)$',
    r'$p \equiv 5\;(8)$',
    r'$p \equiv 17\;(24)$',
    'Case A',
    'NQR mod 7',
    'QR mod 7',
]
values = [cov["3mod4"], cov["5mod8"], cov["17mod24"],
          cov["caseA"], cov["nqr7"], cov["qr7"]]
colors = ['#2f81f7', '#3fb950', '#d29922', '#a371f7', '#f47067', '#f778ba']

wedges, _ = ax1.pie(
    values, labels=None, colors=colors,
    startangle=220, counterclock=False,
    wedgeprops={'linewidth': 1.2, 'edgecolor': '#0d1117'},
    radius=0.95
)

# Labels for big slices only
for i, (w, v) in enumerate(zip(wedges, values)):
    pct = 100 * v / total
    mid = (w.theta2 + w.theta1) / 2
    if pct > 8:
        r = 0.58
        x = r * np.cos(np.radians(mid))
        y = r * np.sin(np.radians(mid))
        ax1.text(x, y, f'{pct:.0f}%', ha='center', va='center',
                fontsize=12, color='white', fontweight='bold')

# Legend below
legend_labels = [f'{labels[i]}  ({100*values[i]/total:.1f}%)'
                 for i in range(len(labels))]
leg = ax1.legend(wedges, legend_labels, loc='upper center',
          bbox_to_anchor=(0.5, -0.02), fontsize=7.5, frameon=False,
          ncol=2, columnspacing=0.8, labelcolor=DIMMED)

ax1.set_title('Prime Classification', pad=8, fontsize=13,
              color=WHITE, fontweight='bold')
ax1.text(0, -1.35, '78,498 primes to $10^6$', ha='center',
         fontsize=8, color=DIMMED, style='italic')

# ============================================================================
# PANEL 2: A-value distribution - bar chart
# ============================================================================
ax2 = fig.add_subplot(gs[1])
ax2.set_facecolor(CARD)

A_dist = {int(k): v for k, v in stats["A_dist"].items()}
total_qr7 = stats["coverage"]["qr7"]
A_vals = sorted(A_dist.keys())
counts = [A_dist[a] for a in A_vals]
pcts = [100 * c / total_qr7 for c in counts]

bars = ax2.bar(range(len(A_vals)), pcts, color=CYAN, edgecolor='#1f6feb',
               linewidth=0.5, alpha=0.85)

# Labels on top of significant bars
for i, (bar, pct) in enumerate(zip(bars, pcts)):
    if pct > 2:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.8,
                f'{pct:.1f}%', ha='center', va='bottom', fontsize=8,
                color=CYAN, fontweight='bold')

# Cumulative line on twin axis
ax2t = ax2.twinx()
cum = np.cumsum(pcts)
ax2t.plot(range(len(A_vals)), cum, 'o-', color=GOLD, markersize=3.5,
          linewidth=1.8, alpha=0.9, zorder=5)
ax2t.set_ylim(0, 115)
ax2t.set_ylabel('Cumulative %', color=GOLD, fontsize=9)
ax2t.tick_params(axis='y', colors=GOLD, labelsize=8)
ax2t.axhline(y=100, color=GOLD, linestyle=':', linewidth=0.7, alpha=0.5)
ax2t.spines['right'].set_color(GOLD)
ax2t.spines['right'].set_alpha(0.3)

ax2.set_xticks(range(len(A_vals)))
ax2.set_xticklabels([str(a) for a in A_vals], rotation=45, ha='right',
                     fontsize=7, color=DIMMED)
ax2.set_ylabel('% of hard primes solved', color=CYAN, fontsize=9)
ax2.tick_params(axis='y', colors=CYAN, labelsize=8)
ax2.set_xlim(-0.6, len(A_vals) - 0.4)

ax2.set_title('Gateway Parameter $A$ Distribution', pad=8,
              fontsize=13, color=WHITE, fontweight='bold')

# Style spines
for spine in ax2.spines.values():
    spine.set_color('#30363d')
for spine in ['top', 'right']:
    ax2.spines[spine].set_visible(False)
ax2t.spines['top'].set_visible(False)
ax2t.spines['left'].set_visible(False)

# ============================================================================
# PANEL 3: Max A stabilisation
# ============================================================================
ax3 = fig.add_subplot(gs[2])
ax3.set_facecolor(CARD)

scale_data = {
    1_000: 7, 2_000: 7, 5_000: 23, 10_000: 31,
    20_000: 47, 50_000: 47, 100_000: 59, 200_000: 59,
    500_000: 79, 1_000_000: 79,
    10_000_000: 167, 100_000_000: 239, 1_000_000_000: 239,
    10_000_000_000: 359, 100_000_000_000: 359,
}
prog = stats.get("max_A_progressive", {})
for k, v in prog.items():
    scale_data[int(k)] = v

limits = sorted(scale_data.keys())
max_As = [scale_data[l] for l in limits]

ax3.semilogx(limits, max_As, 'o-', color=PINK, markersize=6,
             linewidth=2, markeredgecolor='#d63384', markeredgewidth=1,
             zorder=5)

# Earlier plateau (239) and current maximum (359)
ax3.axhline(y=239, color=GOLD, linestyle=':', linewidth=0.9, alpha=0.4)
ax3.text(1.5e5, 245, '$A = 239$', va='bottom', ha='left', fontsize=8,
         color=GOLD, style='italic', alpha=0.7)
ax3.axhline(y=359, color=GOLD, linestyle='--', linewidth=1, alpha=0.6)
ax3.text(2.5e11, 365, '$A = 359$', va='bottom', ha='right', fontsize=9,
         color=GOLD, style='italic')

# Shade the flat 10^10 -> 10^11 region
ax3.axvspan(1e10, 1e11, alpha=0.08, color=PINK)
ax3.text(3.2e10, 25, 'Stable', fontsize=9, color=PINK, alpha=0.7,
         fontweight='bold')

ax3.set_ylabel('Max $A$ needed', color=PINK, fontsize=9)
ax3.set_ylim(0, 400)
ax3.set_xlim(500, 3e11)
ax3.tick_params(axis='y', colors=PINK, labelsize=8)
ax3.tick_params(axis='x', colors=DIMMED, labelsize=8)

xticks = [1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11]
ax3.set_xticks(xticks)
ax3.set_xticklabels(['$10^3$', '', '$10^5$', '', '$10^7$', '',
                      '$10^9$', '', '$10^{11}$'], color=DIMMED)
ax3.grid(True, alpha=0.1, color='#30363d')

ax3.set_title('Max $A$ Stabilisation', pad=8,
              fontsize=13, color=WHITE, fontweight='bold')

for spine in ax3.spines.values():
    spine.set_color('#30363d')
for spine in ['top', 'right']:
    ax3.spines[spine].set_visible(False)

# ============================================================================
# Main title
# ============================================================================
fig.suptitle(
    r'Erdős–Straus Conjecture: Gateway Decompositions to $10^{11}$',
    fontsize=17, color=WHITE, fontweight='bold', y=1.06
)
fig.text(0.5, 0.99,
         '32 values of $A$ resolve all 4,118,054,813 primes  |  max $A$ = 359  |  key insight: $d \\mid N^2$',
         ha='center', fontsize=10, color=GOLD, style='italic')

plt.savefig('banner.png', bbox_inches='tight', facecolor=BG, edgecolor='none', pad_inches=0.3)
print("Banner saved: banner.png")
plt.close()
