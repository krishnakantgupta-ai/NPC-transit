"""Figure 2: leave-one-out predicted versus measured ln(q). Usage: python3 analysis/make_figure2.py results/reg_stokes.properties.loo_predictions.csv figure2.png"""
import csv, sys, math, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
rows = list(csv.DictReader(open(sys.argv[1])))
x = [float(r['y_true']) for r in rows]; y = [float(r['y_loo_pred']) for r in rows]; n = [r['id'] for r in rows]
rm = math.sqrt(sum((a-b)**2 for a, b in zip(x, y))/len(x))
fig, ax = plt.subplots(figsize=(4.8, 4.6), dpi=300); lo, hi = -5.4, 3.4
ax.plot([lo, hi], [lo, hi], color="#888", lw=1, zorder=1)
ax.fill_between([lo, hi], [lo-rm, hi-rm], [lo+rm, hi+rm], color="#d9e6f5", alpha=.7, lw=0, zorder=0, label=f"±1 LOO RMSE ({rm:.2f})")
ax.scatter(x, y, s=28, color="#1f4e79", zorder=3)
for xi, yi, ni in zip(x, y, n): ax.annotate(ni, (xi, yi), textcoords="offset points", xytext=(6, -10), fontsize=7)
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
ax.set_xlabel("Measured ln(q)"); ax.set_ylabel("Leave-one-out predicted ln(q)"); ax.legend(frameon=False, fontsize=7, loc="upper left")
for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout(); plt.savefig(sys.argv[2] if len(sys.argv) > 2 else "figure2.png")
