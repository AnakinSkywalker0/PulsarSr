"""Development diagnostic only; do not use held-out outputs to retune rules."""
from pathlib import Path
import json
root = Path(__file__).resolve().parents[1]
data = json.loads((root / "reference/confirmation_development.json").read_text())
for scenario in ["broadband-rfi", "pulsar", "pulsar-colored", "pulsar-rfi"]:
    print(scenario)
    for row in [r for r in data["observations"] if r["scenario"] == scenario and r["false_accepted"]][:3]:
        print("truth", row["period_s"], row["dm"])
        for c in [c for c in row["diagnostics"] if c["accepted"]]:
            print("candidate", c["profile"]["period_s"], c["refined_dm"],
                  "held", c["validation"]["snr"], "zero", c["zero_dm_validation"]["snr"],
                  "contrast", c["dm_contrast"], "accepted", c["accepted"],
                  "reasons", c["rejection_reasons"])
