"""Plot recovery fractions with Wilson 95% intervals and save a compact summary.
Run with PRESTO's WSL Python (numpy/matplotlib already available).
"""
from pathlib import Path
import json
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parents[1]
data = json.loads((root / "reference/recovery.json").read_text())
rows = data["trials"]
groups = []
for low, high in [(0, 3), (3, 6), (6, 10), (10, 20), (20, 1000)]:
    group = [r for r in rows if r["amplitude"] > 0 and low <= r["expected_matched_filter_snr"] < high]
    if not group:
        continue
    n, k = len(group), sum(r["recovered"] for r in group)
    p, z = k/n, 1.96
    center = (p+z*z/(2*n))/(1+z*z/n)
    half = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    groups.append({"snr_low": low, "snr_high": high, "n": n, "recovered": k,
                   "fraction": p, "ci_low": max(0, center-half), "ci_high": min(1, center+half),
                   "snr_mean": sum(r["expected_matched_filter_snr"] for r in group)/n})
noise = [r for r in rows if r["amplitude"] == 0]
summary = {"bins": groups, "noise_observations": len(noise),
           "noise_with_candidates": sum(r["candidates"] > 0 for r in noise),
           "noise_candidate_total": sum(r["candidates"] for r in noise),
           "protocol": data["protocol"]}
summary["amplitude_slices"] = [
    {"amplitude": a, "n": len(g), "period_recovered": sum(r["period_recovered"] for r in g),
     "strict_recovered": sum(r["recovered"] for r in g)}
    for a in sorted(set(r["amplitude"] for r in rows)) if a > 0
    for g in [[r for r in rows if r["amplitude"] == a]]
]
(root / "reference/recovery_summary.json").write_text(json.dumps(summary, indent=2))
fig, ax = plt.subplots(figsize=(9, 5))
x = [g["snr_mean"] for g in groups]
y = [g["fraction"] for g in groups]
ax.errorbar(x, y, yerr=[[max(0, g["fraction"]-g["ci_low"]) for g in groups],
                       [max(0, g["ci_high"]-g["fraction"]) for g in groups]], marker="o", capsize=4)
for g in groups:
    ax.annotate(f'{g["recovered"]}/{g["n"]}', (g["snr_mean"], g["fraction"]), xytext=(5, -16), textcoords="offset points")
ax.set(xlabel="Analytic injected Gaussian matched-filter SNR proxy", ylabel="Strict period + DM recovery fraction",
       ylim=(-0.05, 1.08), title="Development injection recovery — Wilson 95% intervals")
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(root / "reference/recovery.png", dpi=150)
print(json.dumps(summary, indent=2))
