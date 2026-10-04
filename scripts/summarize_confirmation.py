"""Summarize frozen held-out confirmation without retuning on its outcomes."""
from pathlib import Path
import json
import math

root = Path(__file__).resolve().parents[1]
data = json.loads((root / "reference/confirmation_held-out.json").read_text())
rows = data["observations"]

def interval(k, n):
    z = 1.96
    p = k/n
    center = (p+z*z/(2*n))/(1+z*z/n)
    half = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    return [max(0, center-half), min(1, center+half)]

slices = []
for name in sorted(set(r["scenario"] for r in rows)):
    group = [r for r in rows if r["scenario"] == name]
    recovered = sum(r["recovered"] for r in group)
    slices.append({"scenario": name, "observations": len(group), "recovered": recovered,
                   "recovery_wilson_95": interval(recovered, len(group)),
                   "observations_with_unmatched_accepted": sum(r["false_accepted"] > 0 for r in group),
                   "accepted_total": sum(r["accepted_candidates"] for r in group)})
controls = [r for r in rows if not r["scenario"].startswith("pulsar")]
false_observations = sum(r["accepted_candidates"] > 0 for r in controls)
strength_slices = []
for name in ["pulsar", "pulsar-colored", "pulsar-rfi"]:
    for amplitude in [0.1, 0.2, 0.4, 0.8]:
        g = [r for r in rows if r["scenario"] == name and r["amplitude"] == amplitude]
        strength_slices.append({"scenario": name, "amplitude": amplitude, "n": len(g),
                                "recovered": sum(r["recovered"] for r in g)})
report = json.loads((root / "output_confirm_gbt/report.json").read_text())
accepted = [c for c in report["confirmation"]["candidates"] if c["accepted"]]
result = {"rules": data["rules"], "protocol": data["protocol"], "split": data["split"],
          "slices": slices, "strength_slices": strength_slices,
          "negative_observations": len(controls), "negative_observations_with_accepted": false_observations,
          "false_observation_wilson_95": interval(false_observations, len(controls)),
          "zero_false_exact_one_sided_95_upper": 1-0.05**(1/len(controls)) if false_observations == 0 else None,
          "gbt_accepted": [{"period_ms": c["profile"]["period_s"]*1000, "dm": c["refined_dm"],
                            "validation_snr": c["validation"]["snr"], "dm_contrast": c["dm_contrast"],
                            "time_part_snr": [p["snr"] for p in c["time_parts"]],
                            "subband_snr": [p["snr"] for p in c["subbands"]]} for c in accepted],
          "interpretation": "New random realizations on same template grid; 10 candidates, 32 bins vs GBT CLI 5/64; not general astrophysical false-alarm calibration. Positive unmatched leads can be harmonics, period/DM estimation errors or false positives. Zero negatives observed is not zero population risk. Low-DM/broad pulses may be rejected by conservative vetoes."}
(root / "reference/confirmation_summary.json").write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
