"""Repeat the same Rust workload at 1 and 8 threads; record timing and hardware.

This is a scaling measurement, not an end-to-end PRESTO speed comparison.
"""
from pathlib import Path
import hashlib
import json
import platform
import statistics
import subprocess
import time

root = Path(__file__).resolve().parents[1]
exe = root / "target/release/psr.exe"
dataset = root / "GBT_Lband_PSR.fil"
groups = []
for threads in [1, 8]:
    times = []
    reports = []
    cmd = [str(exe), "search", str(dataset), "--dm-max", "500", "--dm-step", "2",
           "--threads", str(threads), "--top", "2", "--out", str(root / f"output_bench_{threads}")]
    for repeat in range(4):
        start = time.perf_counter()
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        elapsed = time.perf_counter()-start
        report = json.loads((root / f"output_bench_{threads}/report.json").read_text())
        if repeat:  # discard first run as warm-up
            times.append(elapsed)
            reports.append(report["folded"])
    groups.append({"threads": threads, "command": cmd, "wall_seconds": times,
                   "median": statistics.median(times), "stdev": statistics.stdev(times),
                   "folded": reports[0], "repeat_results_identical": all(r == reports[0] for r in reports)})
cpu = subprocess.run(["powershell", "-NoProfile", "-Command",
                      "Get-CimInstance Win32_Processor | Select-Object Name,NumberOfCores,NumberOfLogicalProcessors | ConvertTo-Json"],
                     capture_output=True, text=True, check=True).stdout
result = {"platform": platform.platform(), "cpu": json.loads(cpu),
          "rustc": subprocess.run(["rustc", "--version"], capture_output=True, text=True, check=True).stdout.strip(),
          "dataset_sha256": hashlib.sha256(dataset.read_bytes()).hexdigest(),
          "protocol": "Release; 251 DMs step 2; default RFI; top 2; one warm-up then three timed repeats; warm file cache; includes plots; same channel summation order at all thread counts",
          "thread_results_identical": groups[0]["folded"] == groups[1]["folded"], "groups": groups}
(root / "reference/benchmark.json").write_text(json.dumps(result, indent=2))
print(json.dumps({"cpu": result["cpu"], "timings": [{k: g[k] for k in ["threads", "median", "stdev"]} for g in groups],
                  "thread_results_identical": result["thread_results_identical"]}, indent=2))
