"""Generate shareable plots and a factual evidence summary from saved runs.
Run with the local PRESTO environment's Python (matplotlib/numpy).
"""
from pathlib import Path
import json
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parents[1]
out = root / "docs/assets"
out.mkdir(parents=True, exist_ok=True)
report = json.loads((root / "output_gbt/report.json").read_text())
# Identify the tutorial target for this comparison only. This does not alter
# search ranking or claim the largest statistic necessarily belongs to a pulsar.
matches = [f for f in report["folded"] if abs(f["profile"]["period_s"]-0.004621638) < 0.00001]
if not matches:
    raise RuntimeError("Tutorial target absent from folded candidates")
best = max(matches, key=lambda f: f["profile"]["snr"])
profile = best["profile"]
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(np.arange(len(profile["bins"]))/len(profile["bins"]), profile["bins"])
ax.set(xlabel="Pulse phase", ylabel="Mean summed power",
       title=f'GBT tutorial pulsar: {profile["period_s"]*1000:.6f} ms, DM {best["refined_dm"]:.3f}')
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(out / "pulse_profile.png", dpi=150)
plt.close(fig)

curve = best["dm_curve"]
fig, ax = plt.subplots(figsize=(9, 4))
ax.plot([x[0] for x in curve], [x[1] for x in curve], marker=".")
ax.set(xlabel="DM (pc cm^-3)", ylabel="Uncalibrated profile SNR", title="Local DM response")
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(out / "dm_response.png", dpi=150)
plt.close(fig)

presto_lines = (root / "reference/presto/gbt_fold_4.62ms_Cand.pfd.bestprof").read_text().splitlines()
presto_profile = np.array([float(line.split()[1]) for line in presto_lines if line and not line.startswith("#")])
rust_profile = np.array(profile["bins"])
# Normalize shape and align phase by the best circular integer-bin shift.
presto_norm = (presto_profile-presto_profile.mean())/presto_profile.std()
rust_norm = (rust_profile-rust_profile.mean())/rust_profile.std()
correlations = [float(np.corrcoef(rust_norm, np.roll(presto_norm, shift))[0, 1]) for shift in range(len(rust_norm))]
shift = int(np.argmax(correlations))
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(rust_norm, label="Rust")
ax.plot(np.roll(presto_norm, shift), label=f"PRESTO (phase shift {shift} bins)", alpha=0.7)
ax.set(xlabel="Phase bin", ylabel="Standardized profile", title="Profile-shape comparison (different folding algorithms/durations)")
ax.legend()
fig.tight_layout()
fig.savefig(out / "presto_comparison.png", dpi=150)
plt.close(fig)
presto_period = float(re.search(r"# P_topo \(ms\)\s*=\s*([\d.]+)", "\n".join(presto_lines)).group(1))
presto_dm = float(re.search(r"# Best DM\s*=\s*([\d.]+)", "\n".join(presto_lines)).group(1))
evidence = {
    "rust_period_ms": profile["period_s"]*1000, "rust_dm": best["refined_dm"],
    "presto_period_ms": presto_period, "presto_dm": presto_dm,
    "period_difference_microseconds": (profile["period_s"]*1000-presto_period)*1000,
    "profile_shape_correlation_after_phase_alignment": max(correlations),
    "rust_rfi_masked_fraction": report["rfi"]["masked_fraction"],
    "search_seconds_single_run": report["elapsed_seconds"],
    "dm_trials": report["trial_count"], "dm_step": report["dm_step"],
    "reference": json.loads((root / "reference/presto/comparison.json").read_text()),
    "limits": "PRESTO fold uses 460800 samples/subband folding; Rust uses shared search length/integer channel delays. Scores differ and are not comparable sigmas. Catalogue cross-check and matched-workload speed comparison pending.",
}
(root / "reference/evidence.json").write_text(json.dumps(evidence, indent=2))
print(json.dumps(evidence, indent=2))
