"""Run independent topocentric PRESTO checks under the local WSL environment.

Run: wsl -d Ubuntu -- reference/toolchain/presto/bin/python scripts/reference_run.py
Rust's `psr dedisperse ... --dm 62 --out output_dm62` must be run first.
"""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "reference/toolchain/presto/bin"
OUT = ROOT / "reference/presto"
OUT.mkdir(parents=True, exist_ok=True)
env = dict(os.environ, PATH=str(BIN) + os.pathsep + os.environ.get("PATH", ""))
data = ROOT / "GBT_Lband_PSR.fil"
runs = []

def run(name, *args):
    cmd = [str(BIN / args[0]), *map(str, args[1:])]
    start = time.perf_counter()
    with (OUT / f"{name}.log").open("w") as log:
        result = subprocess.run(cmd, cwd=OUT, env=env, stdout=log, stderr=subprocess.STDOUT)
    runs.append({"command": cmd, "seconds": time.perf_counter()-start, "returncode": result.returncode})
    (OUT / "commands.json").write_text(json.dumps(runs, indent=2))
    if result.returncode:
        raise RuntimeError(f"{name} failed; see {OUT / (name + '.log')}")

run("header", "readfile", data)
run("rfi", "rfifind", "-time", "1", "-o", "gbt_rfi", data)
run("dm62", "prepdata", "-nobary", "-noclip", "-dm", "62", "-numout", "530000", "-o", "gbt_dm62", data)
run("fft", "realfft", "gbt_dm62.dat")
run("candidates", "accelsearch", "-zmax", "0", "-flo", "15", "-numharm", "16", "gbt_dm62.fft")
run("fold", "prepfold", "-topo", "-noclip", "-noxwin", "-p", "0.004621638", "-dm", "62", "-n", "64", "-nsub", "32", "-npart", "32", "-nopdsearch", "-o", "gbt_fold", data)

rust = np.fromfile(ROOT / "output_dm62/series.f32", dtype="<f4")
presto = np.fromfile(OUT / "gbt_dm62.dat", dtype="<f4")
n = min(len(rust), len(presto))
# Exclude end padding/buffer boundaries; no fitted time shift in the headline check.
a, b = rust[:n-2400].astype(float), presto[:n-2400].astype(float)
correlation = float(np.corrcoef(a, b)[0, 1])
slope, intercept = np.polyfit(a, b, 1)
residual_rms = float(np.sqrt(np.mean((b-(slope*a+intercept))**2)))
comparison = {
    "dataset_sha256": hashlib.sha256(data.read_bytes()).hexdigest(),
    "presto_version": "6.1.0 (conda-forge presto-pulsar)",
    "comparison": "DM62, nobary, noclip, unmasked; aligned at start, final 2400 samples excluded",
    "samples_compared": len(a), "correlation": correlation,
    "presto_equals_rust_scale": float(slope), "offset": float(intercept),
    "residual_rms": residual_rms, "passed_correlation_gate": correlation > 0.99,
}
(OUT / "comparison.json").write_text(json.dumps(comparison, indent=2))
print(json.dumps(comparison, indent=2))
