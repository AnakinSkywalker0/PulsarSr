# Building PulsarSr: correctness before speed

## Problem and baseline

Radio pulses arrive later at lower frequencies. The first challenge is recovering a periodic signal after aligning channels at an unknown dispersion measure. PulsarSr begins with a brute-force numerical implementation rather than an optimized approximation; PRESTO supplies an independent baseline.

The small GBT tutorial file is particularly useful: its packed 4-bit samples, descending channel order and known millisecond signal test the entire path. Rust's raw DM62 time series agrees exactly with PRESTO in the recorded comparable interval.

## Design

Seven crates isolate file format, cleaning, dedispersion, search, folding, synthesis and CLI. Channel-major numerical storage makes delayed channel slices contiguous. Each DM job has independent FFT buffers and preserves channel summation order; parallelization changes scheduling rather than reductions. JSON reports and synthetic truth sidecars make results inspectable.

## Bugs found through testing

### Spectral normalization on noise-free data

A pure-sine fixture initially ranked tiny FFT leakage at high frequencies above the real peak. Local medians were effectively zero, amplifying numerical residues. A floor tied to the global average spectrum prevents dividing by a nearly zero noise estimate. The pure-sine regression now passes.

### Off-bin harmonic sums

The first development recovery curve plateaued at 15/40 at the strongest amplitude. Inspection showed true signals represented at half their periods: direct `k × harmonic` lookup drifts away from higher harmonics when the fundamental is between FFT bins. Zero padding and local peak interpolation before harmonic sampling recover the intended fundamental in the strong injection slice. A seeded broadband-noise off-bin regression guards the change.

The revised development results reuse the same seeds and therefore are diagnostic evidence, not an independent test-set improvement claim.

### Profile refinement objective

Maximizing peak/off-pulse RMS alone favored unstable off-pulse variance in local period/DM grids. A flat-profile reduced chi-square objective now drives refinement; peak SNR remains an explicitly uncalibrated visualization/ranking measure. Spectral interference can still have high chi-square, so neither statistic alone is astrophysical proof.

### Burst RFI masking floor

Severe-RFI recovery was 0/32. The cause was in the mask, not the confirmation rules: the robust scale was floored at 1% of the median, about 1.0 on data near 100, larger than the true scatter of block means (about 0.0125). A broadband burst shifting a block mean by about 6 stayed under the 8-sigma threshold and was never masked. After dedispersion it produced strong 15-20 Hz candidates near DM 80 that crowded out the pulsar and caused false accepts. Lowering the floor to 0.01% raised recovery to 20/32 on fresh seeds (base 60000), matching the clean slice, with negatives still 0/160 and the GBT result unchanged. Many wrong-harmonic misses seen before the fix were a symptom of the same burst.

Remaining misses were then classified: signals too weak to detect (recovery by injected amplitude 0.1/0.2/0.4/0.8 was 1/8, 5/8, 7/8, 8/8 clean), and correctly found low-DM or broad-pulse pulsars vetoed by the zero-DM contrast and dispersion-width rules. Relaxing those rules trades directly against false-alarm control, so it is left as an explicit decision for a fresh held-out set.
## Interpretation

PRESTO folds the observation at ~4.621376 ms topocentrically, slightly different from the tutorial's quoted 4.621638 ms. Comparing values requires keeping timing conventions explicit. Rust does not barycenter, and its DM/profile refinement differs from PRESTO's subband method.

Strong periodic interference is present. The current pipeline returns candidate leads and plots; a robust calibrated confirmation policy remains unfinished. Noise-only controls make this limitation visible rather than hiding it behind injection recovery.

## Next work

1. Independently check the PRESTO-suggested catalogue identification and timing parameters.
2. Split-data confirmation and fresh held-out evaluation are implemented in [confirmation.md](confirmation.md). Improve the documented low-DM/harmonic/RFI failures using development data, then reserve another held-out seed set.
3. Add realistic RFI, colored noise, fractional-sample and short-period stress slices.
4. Compare full candidate lists and injected recovery with PRESTO under matched settings.
5. Profile and implement subband/tree dedispersion, checking against brute force.
6. Run matched PRESTO throughput measurements before claiming a speed advantage.
7. Choose optional PSRFITS/single-pulse/acceleration extensions after the core validation is stronger.
