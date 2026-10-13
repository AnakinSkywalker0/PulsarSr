# Candidate confirmation — split-confirm-v1

## Changes

The search now preserves spectral detections across DM trials and harmonic-family members, generates full-range DM response diagnostics, and distinguishes **accepted follow-up leads** from rejected candidates. `psr fold` is available for a specified period/DM without searching or declaring a detection.

Report schema 2 records the input SHA-256, available Git revision (null before a repository exists), split sizes, exact rules, sample counts, on/off statistics, rejection reasons, four held-out time profiles and four frequency-subband profiles. Discovery and validation reports are preserved; rejected candidates are not silently discarded.

## Parameter selection and independence

1. Reserve the first usable half for FFT candidate generation and period/DM/window refinement.
2. Insert a maximum-channel-delay guard gap so shifted discovery samples do not overlap raw-channel samples used by validation at another DM.
3. Select the phase and pulse width using discovery data only, with power-of-two windows no wider than one quarter of a rotation.
4. Measure the selected window in the second half without changing candidate period, DM, phase or width.
5. Apply the same fixed window and absolute phase to four validation time segments and four validation subbands.

RFI masking remains an unsupervised full-observation preprocessing step, so strict independence is conditional on that mask. Colored/correlated noise and mask selection can invalidate nominal independent-sample standard errors. Empirical controls, rather than the SNR number alone, are therefore required.

## Frozen acceptance criteria

All must pass:

- Held-out on/off pulse contrast SNR ≥ **6**.
- At least **3/4** time parts with the same-window SNR ≥ **1**.
- At least **3/4** subbands with same-window SNR ≥ **1**.
- Nonzero refined DM and validation SNR at least **1.25×** the zero-DM control SNR (denominator floored at 1).
- Band delay spans at least **1.5 selected pulse widths**.
- No stronger rationally related zero-DM discovery feature, considering ratios of integers 1–16.

The zero-DM rejection control is conservatively optimized in period/phase/width on validation data. This can only make acceptance harder; its statistic is not an independent significance measurement. The rational-frequency veto can reject coincident astrophysical signals as well as RFI. Broad or low-DM sources can fail the dispersion-width check. These choices are documented sensitivity trade-offs, not universal pulsar criteria.

SNR here is the mean on/off sample difference divided by the estimated standard error of that difference. It is not PRESTO sigma or a globally calibrated false-alarm probability.

## DM diagnostics

Each trial contains its actual discovery spectral score/frequency if detected above threshold; missing detections are **null**, not zero. All harmonic members remain in JSON.

The full-range folded response uses first-half per-channel phase histograms shifted/interpolated by delay. It is a visualization approximation: phase-bin interpolation and edge truncation differ from direct dedispersion. Acceptance and refined DM use directly dedispersed discovery/validation samples, not that approximation.

## Held-out experiment

Rules were revised on development seeds based at 20000, frozen, then run on new seeds based at 40000. There are 256 held-out observations, 32 per scenario. The pulse template grid repeats development periods/DMs; the random realizations are new, but this is not out-of-distribution generalization.

| Scenario | Strict recoveries | Observations with unmatched accepted leads |
| --- | --- | --- |
| White noise | — | 0/32 |
| Colored noise | — | 0/32 |
| Broadband zero-DM periodic RFI | — | 0/32 |
| Narrowband sinusoidal RFI | — | 0/32 |
| Burst plus narrowband RFI | — | 0/32 |
| Pulsar + white noise | 17/32 | 2/32 |
| Pulsar + colored noise | 16/32 | 6/32 |
| Pulsar + burst/narrowband RFI | 0/32 | 7/32 |

**No accepted candidates occurred in 160 negative observations.** The exact one-sided 95% upper bound is about **1.86% per observation** under this control mixture; each 32-observation subgroup has much larger uncertainty. Do not interpret the observation of zero as a zero-risk claim.

Positive unmatched leads include wrong harmonics, DM errors and possible noise candidates; they are not all proven noise false alarms. Strict recovery requires fundamental frequency within 1.5 discovery Fourier bins and DM within 6, without harmonic credit.

Strengths vary over amplitudes 0.1/0.2/0.4/0.8 and duty cycles 0.03/0.05/0.1. These are aggregates over strength. Per-strength counts and Wilson intervals are in `reference/confirmation_summary.json`.

The colored stress generator adds an AR(1)-filtered copy of the original channel fluctuations (which include the pulse for positive cases). This is a controlled distortion/noise slice, not a physical radio-noise simulator. Burst RFI can fill the limited discovery budget with aliases and spoil DM estimates. Its complete strict-recovery failure is an important unresolved limitation.

## GBT observation

The new split analysis accepts one of five followed-up candidates:

- **P ≈ 4.621398 ms**, **DM ≈ 62.389**
- Held-out fixed-window SNR **≈ 34.18**
- The strongest ~34.77 ms interference feature is rejected.

Those values are selected from less discovery data and are not replacements for the earlier full-data baseline measurements (P 4.621383 ms, DM 62.014). Historical PRESTO comparisons remain recorded separately.

## Reproduce

```text
cargo run --release -p synth --example confirmation -- development reference/confirmation_development.json
cargo run --release -p synth --example confirmation -- held-out reference/confirmation_held-out.json
target\release\psr.exe search GBT_Lband_PSR.fil --dm-max 500 --threads 8 --out output_confirm_gbt
python scripts\summarize_confirmation.py
target\release\psr.exe fold GBT_Lband_PSR.fil --period 0.004621398 --dm 62.389 --out output_fold
```

The synthetic evaluation uses **10** families and **32** bins; the displayed GBT command uses default **5** families and **64** bins. Results depend on these budgets. Do not compare them as identical statistical workloads.

## Still required

Realistic RFI mitigation, better rational/harmonic family recovery, robust correlated-noise errors, more independent negative observations, a new held-out set after future tuning, ATNF/epoch checks and complete PRESTO sensitivity comparisons. Acceptance currently prioritizes a narrower set of persistent dispersed pulses; it is not astrophysical proof.

## Addendum: burst-RFI mask fix

Root cause of the 0/32 severe-RFI recovery: `rfi::mask` floored its robust scale at 1% of the median, about 1.0 for data near 100, so a +6 block-mean burst fell under the 8-sigma threshold and was never masked. The dedispersed burst then produced low-frequency (15-20 Hz) candidates at DM ~80 that crowded out the pulsar. Floor reduced to 0.01% with a regression test.

Results (strict recovery / observations with unmatched accepted leads), development seeds 20000 then fresh `held-out-2` seeds 60000: pulsar-rfi 0/32 -> 20/32 (1/32 false) on development and 20/32 (1/32 false) on held-out-2; negatives 0/160 on held-out-2; GBT unchanged. Remaining misses are harmonic/DM errors shared with the clean slice. Raw output: `reference/confirmation_held-out-2.json`.
