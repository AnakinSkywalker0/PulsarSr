"""Bootstrap a workspace-local Linux PRESTO environment (run via WSL).

Uses the upstream recommended conda-forge presto-pulsar distribution.
No administrator rights or global shell configuration changes are required.
"""
from pathlib import Path
import subprocess
import tarfile
import urllib.request

root = Path(__file__).resolve().parents[1]
tools = root / "reference" / "toolchain"
tools.mkdir(parents=True, exist_ok=True)
archive = tools / "micromamba.tar.bz2"
binary = tools / "bin" / "micromamba"
if not binary.exists():
    urllib.request.urlretrieve("https://micro.mamba.pm/api/micromamba/linux-64/latest", archive)
    with tarfile.open(archive) as tar:
        tar.extract("bin/micromamba", tools, filter="data")
    binary.chmod(0o755)
subprocess.run([
    str(binary), "create", "-y", "--root-prefix", str(tools / "mamba"),
    "--prefix", str(tools / "presto"), "-c", "conda-forge", "presto-pulsar",
], check=True)
