"""Exploratory paired comparison of two leave-one-out regressions, plus leave-one-protein-out sensitivity.
Usage: python3 analysis/paired_comparison.py results/reg_stokes.properties.loo_predictions.csv results/reg_mass.properties.loo_predictions.csv
Reports, per protein, |error of second model| - |error of first model| (positive = first model better)."""
import csv, random, statistics as st, sys
def load(f): return {r['id']: (float(r['y_true']), float(r['y_loo_pred'])) for r in csv.DictReader(open(f))}
a, b = load(sys.argv[1]), load(sys.argv[2])
ids = sorted(a); d = [abs(b[i][1]-b[i][0]) - abs(a[i][1]-a[i][0]) for i in ids]
for i, x in zip(ids, d): print(f"{i:12s} {x:+.3f}")
print("mean difference", round(st.mean(d), 3), "| first model better in", sum(x > 0 for x in d), "of", len(d))
random.seed(1); bs = sorted(st.mean(random.choices(d, k=len(d))) for _ in range(5000))
print("95% bootstrap CI of mean difference:", round(bs[125], 3), round(bs[4874], 3))
print("leave-one-protein-out mean differences:", [round(st.mean(d[:k]+d[k+1:]), 3) for k in range(len(d))])
