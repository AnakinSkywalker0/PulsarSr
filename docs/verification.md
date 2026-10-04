# Verification record

**Update:** split-data confirmation, fresh held-out controls and independent reader fixtures are documented in [confirmation.md](confirmation.md). Historical full-data comparisons below retain their original settings and output files.

## Dataset

`GBT_Lband_PSR.fil`: 25,488,365 bytes, SHA-256
`0ef6a018fc717bb21af368360a4c90f8e90cf9cc08705b728836c150047205e8`.

Rust and PRESTO agree on 365-byte header, 531,000 samples, 96 channels, 4-bit data, 72 µs sample time, 1447.5 MHz first channel, −1 MHz channel offset, MJD 53010.48482638889. PRESTO additionally derives telescope/sky/date metadata; raw keyword JSON is not literally the same display as `readfile`.

## PRESTO independent run

PRESTO 6.1.0 from conda-forge `presto-pulsar`, installed using workspace-local micromamba in WSL Ubuntu. No Docker daemon was needed.

From the Windows workspace:

```text
wsl -d Ubuntu -- python3 /mnt/d/PulsarSr/scripts/setup_presto.py
target\release\psr.exe dedisperse GBT_Lband_PSR.fil --dm 62 --out output_dm62
wsl -d Ubuntu -- /mnt/d/PulsarSr/reference/toolchain/presto/bin/python /mnt/d/PulsarSr/scripts/reference_run.py
wsl -d Ubuntu -- /mnt/d/PulsarSr/reference/toolchain/presto/bin/python /mnt/d/PulsarSr/scripts/summarize.py
```

The script saves every actual command and duration in `reference/presto/commands.json`, plus logs and binary outputs. It runs `readfile`, `rfifind`, `prepdata`, `realfft`, `accelsearch` and `prepfold`.

### Time series

`prepdata -nobary -noclip -dm 62 -numout 530000`, no mask. Rust similarly sums unmasked channels. Compare aligned starts, excluding final 2,400 samples to avoid buffer/padding boundaries. Over 527,600 samples:

- Pearson correlation: **1.0**
- Fitted PRESTO/Rust scale: **0.9999999999999996**
- Residual RMS after fit: **3.41 × 10⁻¹³** (fit roundoff)

This independently exercises real 4-bit packed decoding, channel order, physical units and integer delay conventions. See `reference/presto/comparison.json`.

### Search and folding

PRESTO `accelsearch -zmax 0 -flo 15 -numharm 16` on DM62 ranks:

1. 4.62139(3) ms / 216.385(2) Hz, notes `PSR J1643-1224?`
2. ~34.77 ms / 28.76 Hz interference-like feature

Rust detects both but ranks the interference first under its different statistic. A candidate's spectral rank is not its astrophysical confirmation.

PRESTO raw-data topocentric fold: period **4.62137647351207 ms**, DM **62.318**. Rust's final full-range search measurements and phase-aligned profile correlation are in `reference/evidence.json`. The catalogue association suggested by PRESTO has not yet been independently checked against ATNF.

Differences that affect agreement:

- PRESTO folds 460,800 samples into subbands/subintegrations; Rust folds the common valid search interval, at integer channel delays.
- PRESTO refines a broader period/DM grid; Rust uses a local frequency grid and ±5 original DM steps at half-step resolution.
- PRESTO's profile statistics and Rust's peak SNR/reduced flat-profile chi-square have different definitions.
- PRESTO's `rfifind` bad fraction is **3.632%** of analyzed channel/time intervals, separately **2.703% padded**. Rust masks **2.484%** of actual samples. Different intervals, denominators and algorithms prevent interpreting those percentages as mask equivalence.

## Injection protocol

`cargo run --release -p synth --example recovery -- reference/recovery.json`

240 paired seeded observations: 40 noise-only controls and 200 Gaussian pulse injections. Each has 16 channels, 16,384 samples, 250 µs sampling. Periods: 23.71/32/47.19/61.37 ms; DMs: 18/42/62/96/118; duty cycles: 0.03/0.05/0.1; amplitudes 0.025/0.05/0.1/0.2/0.4 against unit-variance channel noise. Seeds 1000–1039 repeat across amplitudes for paired comparison.

Search DM 0–120 step 2; frequency 10–1000 Hz; spectral threshold 8; at most 10 sifted candidates. Strict recovery requires the **fundamental**, within 1.5 unpadded Fourier bins and DM within 4. No harmonic credit. Per-trial truth, errors, scores and candidate counts are saved.

Injected SNR proxy is `A sqrt(Nchannel Ntime (duty/2.355) sqrt(pi))`, Gaussian matched-filter energy in white noise. It is not the FFT score or measured folded SNR.

Generate summary and Wilson intervals:

```text
wsl -d Ubuntu -- /mnt/d/PulsarSr/reference/toolchain/presto/bin/python /mnt/d/PulsarSr/scripts/plot_recovery.py
```

This is **development evaluation**: the data were reused while fixing off-bin harmonic summation. It does not support a held-out generalization claim. Above SNR proxy 20, strict recovery is 50/53; all 40 pure-noise observations have candidates at this permissive threshold. Final confirmation thresholds and false-alarm calibration remain to be designed.

## Checks run

- `cargo test --workspace`: now 14 unit/integration tests passed, including seeded end-to-end confirmation, zero-DM/noise controls and independent literal bit-depth payloads.
- `cargo clippy --workspace --all-targets -- -D warnings`: passed.
- `cargo fmt --all -- --check`: passed after formatting.
- Real-data header, inspect, dedispersion, search and folding commands executed.
- Thread-scaling benchmark repeated with matching outputs at 1/8 threads.

CI configuration exists; remote CI and exhaustive scientific validation are pending.
