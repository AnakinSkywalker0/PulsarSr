//! Discovery/confirmation split: parameter/window selection uses the first half
//! only. The second half is measured without optimizing phase, period, DM or width.
//! Scores assume independent finite-variance samples; colored noise/RFI require
//! empirical controls. An accepted candidate is a follow-up lead, not proof of origin.
use crate::{Candidate, Config};
use serde::Serialize;

#[derive(Clone, Copy)]
pub struct Observation<'a> {
    pub data: &'a [f32],
    pub nsamples: usize,
    pub frequencies: &'a [f64],
    pub tsamp: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct Rules {
    pub version: String,
    pub minimum_validation_snr: f64,
    pub minimum_part_snr: f64,
    pub minimum_time_parts: usize,
    pub minimum_subbands: usize,
    pub minimum_dm_contrast: f64,
    pub minimum_dispersion_widths: f64,
    pub maximum_zero_dm_spectral_ratio: f64,
}
impl Default for Rules {
    fn default() -> Self {
        Self {
            version: "split-confirm-v1".into(),
            minimum_validation_snr: 6.,
            minimum_part_snr: 1.,
            minimum_time_parts: 3,
            minimum_subbands: 3,
            minimum_dm_contrast: 1.25,
            minimum_dispersion_widths: 1.5,
            maximum_zero_dm_spectral_ratio: 1.,
        }
    }
}
#[derive(Clone, Debug, Serialize)]
pub struct Window {
    pub start_bin: usize,
    pub width_bins: usize,
    pub bins: usize,
}
impl Window {
    fn contains(&self, bin: usize) -> bool {
        (bin + self.bins - self.start_bin) % self.bins < self.width_bins
    }
}
#[derive(Clone, Debug, Serialize)]
pub struct Measurement {
    pub snr: f64,
    pub on_samples: usize,
    pub off_samples: usize,
    pub mean_difference: f64,
    pub standard_error: f64,
}
#[derive(Clone, Debug, Serialize)]
pub struct DmPoint {
    pub dm: f64,
    /// No detection above the generation threshold is represented by null, not zero.
    pub spectral_score: Option<f64>,
    pub spectral_frequency_hz: Option<f64>,
    /// Approximation from shifted first-half channel phase histograms. Not used
    /// for parameter refinement or the final acceptance decision.
    pub discovery_profile_response: f64,
}
#[derive(Clone, Debug, Serialize)]
pub struct ConfirmedCandidate {
    pub candidate: Candidate,
    pub refined_dm: f64,
    pub profile: fold::Profile,
    pub window: Window,
    pub validation: Measurement,
    pub zero_dm_validation: Measurement,
    pub zero_dm_period_s: f64,
    pub zero_dm_window: Window,
    pub dm_contrast: f64,
    pub dispersion_widths: f64,
    pub zero_dm_related_spectral_ratio: f64,
    pub time_parts: Vec<Measurement>,
    pub subbands: Vec<Measurement>,
    pub time_profiles: Vec<Vec<f32>>,
    pub subband_profiles: Vec<Vec<f32>>,
    pub dm_response: Vec<DmPoint>,
    pub harmonic_members: Vec<Candidate>,
    pub accepted: bool,
    pub rejection_reasons: Vec<String>,
}
#[derive(Clone, Debug, Serialize)]
pub struct Outcome {
    pub rules: Rules,
    pub discovery_samples: usize,
    pub validation_start_sample: usize,
    pub validation_samples: usize,
    pub common_samples: usize,
    pub trial_count: usize,
    pub candidates: Vec<ConfirmedCandidate>,
}

fn phase_bin(t: usize, tsamp: f64, period: f64, bins: usize) -> usize {
    (((t as f64 * tsamp / period).fract() * bins as f64) as usize).min(bins - 1)
}
/// A fixed on/off window measurement; no phase or width optimization.
pub fn measure(x: &[f32], offset: usize, tsamp: f64, period: f64, window: &Window) -> Measurement {
    let mut n = [0usize; 2];
    let mut sum = [0.; 2];
    let mut sum2 = [0.; 2];
    for (i, &v) in x.iter().enumerate() {
        let k = usize::from(window.contains(phase_bin(i + offset, tsamp, period, window.bins)));
        n[k] += 1;
        sum[k] += v as f64;
        sum2[k] += (v as f64).powi(2);
    }
    if n.iter().any(|&v| v < 2) {
        return Measurement {
            snr: 0.,
            on_samples: n[1],
            off_samples: n[0],
            mean_difference: 0.,
            standard_error: 0.,
        };
    }
    let difference = sum[1] / n[1] as f64 - sum[0] / n[0] as f64;
    let variance =
        |k: usize| ((sum2[k] - sum[k] * sum[k] / n[k] as f64) / (n[k] - 1) as f64).max(0.);
    let error = (variance(1) / n[1] as f64 + variance(0) / n[0] as f64).sqrt();
    Measurement {
        snr: if error > 1e-10 {
            difference / error
        } else {
            0.
        },
        on_samples: n[1],
        off_samples: n[0],
        mean_difference: difference,
        standard_error: error,
    }
}
fn select_window(profile: &fold::Profile) -> Window {
    let bins = profile.bins.len();
    let mut best = (
        f64::NEG_INFINITY,
        Window {
            start_bin: 0,
            width_bins: 1,
            bins,
        },
    );
    for width in (0..usize::BITS)
        .map(|i| 1usize << i)
        .take_while(|&w| w <= bins / 4)
    {
        for start in 0..bins {
            let on = (0..width)
                .map(|j| profile.bins[(start + j) % bins] as f64)
                .sum::<f64>();
            let total = profile.bins.iter().map(|&v| v as f64).sum::<f64>();
            let contrast = on / width as f64 - (total - on) / (bins - width) as f64;
            let score = contrast / (1. / width as f64 + 1. / (bins - width) as f64).sqrt();
            if score > best.0 {
                best = (
                    score,
                    Window {
                        start_bin: start,
                        width_bins: width,
                        bins,
                    },
                );
            }
        }
    }
    best.1
}
fn group(
    candidates: &[Candidate],
    resolution: f64,
    limit: usize,
) -> Vec<(Candidate, Vec<Candidate>)> {
    let mut sorted = candidates.to_vec();
    sorted.sort_by(|a, b| b.score.total_cmp(&a.score));
    let mut families: Vec<(Candidate, Vec<Candidate>)> = vec![];
    for c in sorted {
        let related = |a: &Candidate| {
            (1..=16).any(|h| {
                (a.frequency_hz - c.frequency_hz * h as f64).abs() <= resolution
                    || (c.frequency_hz - a.frequency_hz * h as f64).abs() <= resolution
            })
        };
        if let Some((_, members)) = families.iter_mut().find(|(a, _)| related(a)) {
            members.push(c);
        } else if families.len() < limit {
            families.push((c.clone(), vec![c]));
        }
    }
    families
}

/// Analyze already-generated FIRST-HALF candidates with a common delay-trimmed
/// observation. Keeping generation external permits parallel DM FFTs in the CLI.
pub fn analyze(
    observation: Observation<'_>,
    trials: &[f64],
    generated: &[Candidate],
    limit: usize,
    bins: usize,
    rules: &Rules,
) -> Result<Outcome, dedisperse::Error> {
    let Observation {
        data,
        nsamples,
        frequencies,
        tsamp,
    } = observation;
    if trials.is_empty()
        || limit == 0
        || bins < 8
        || frequencies.len() < 4
        || nsamples.checked_mul(frequencies.len()) != Some(data.len())
    {
        return Err(dedisperse::Error(
            "confirmation requires >=4 channels, >=8 bins and valid array/trials".into(),
        ));
    }
    if trials.windows(2).any(|w| w[1] <= w[0]) {
        return Err(dedisperse::Error("DM trials must increase".into()));
    }
    let max_dm = *trials.last().unwrap();
    let delays = dedisperse::delays(frequencies, tsamp, max_dm)?;
    let length = nsamples
        .checked_sub(*delays.iter().max().unwrap())
        .filter(|&n| n >= 512)
        .ok_or_else(|| dedisperse::Error("observation too short for confirmation".into()))?;
    // Gap prevents high-DM discovery shifts reading samples subsequently used
    // for low-DM validation on the same channel.
    let guard = *delays.iter().max().unwrap();
    let discovery = length
        .checked_sub(guard)
        .filter(|&n| n >= 256)
        .ok_or_else(|| dedisperse::Error("DM delay leaves no independent split".into()))?
        / 2;
    let validation_start = discovery + guard;
    let validation = length - validation_start;
    let resolution = 1. / (discovery as f64 * tsamp);
    let step = if trials.len() > 1 {
        trials[1] - trials[0]
    } else {
        1.
    };
    let zero_series =
        dedisperse::apply(data, nsamples, &vec![0; frequencies.len()], Some(discovery))?;
    let zero_candidates = crate::periodicity(
        &zero_series,
        tsamp,
        0.,
        &Config {
            min_hz: 1.,
            max_hz: 0.5 / tsamp,
            threshold: 8.,
            max_candidates: 64,
            birdies: vec![],
        },
    );
    let mut results = vec![];
    for (candidate, members) in group(generated, resolution * 1.5, limit) {
        // Harmonic sums can produce rational-frequency aliases of a local emitter.
        let zero_ratio = zero_candidates
            .iter()
            .filter(|z| {
                (1..=16).any(|m| {
                    (1..=16).any(|n| {
                        (candidate.frequency_hz - z.frequency_hz * m as f64 / n as f64).abs()
                            <= resolution * 1.5
                    })
                })
            })
            .map(|z| z.score / candidate.score.max(1.))
            .fold(0., f64::max);
        let initial = dedisperse::apply(
            data,
            nsamples,
            &dedisperse::delays(frequencies, tsamp, candidate.dm)?,
            Some(length),
        )?;
        let mut best = fold::refine(
            &initial[..discovery],
            tsamp,
            candidate.frequency_hz,
            resolution,
            40,
            bins,
        );
        let mut best_dm = candidate.dm;
        let period = best.period_s;
        for i in -10..=10 {
            let dm = (candidate.dm + i as f64 * step / 2.).clamp(trials[0], max_dm);
            let x = dedisperse::apply(
                data,
                nsamples,
                &dedisperse::delays(frequencies, tsamp, dm)?,
                Some(discovery),
            )?;
            let p = fold::profile(&x, tsamp, period, bins);
            if p.reduced_chi_square > best.reduced_chi_square {
                best = p;
                best_dm = dm;
            }
        }
        let window = select_window(&best);
        let d = dedisperse::delays(frequencies, tsamp, best_dm)?;
        let series = dedisperse::apply(data, nsamples, &d, Some(length))?;
        let validation_measure = measure(
            &series[validation_start..],
            validation_start,
            tsamp,
            period,
            &window,
        );
        let zero = dedisperse::apply(data, nsamples, &vec![0; frequencies.len()], Some(length))?;
        // Compare zero-DM at every phase with the SAME fixed width, conservatively
        // favoring interference. No zero-DM phase adjustment boosts the accepted score.
        // Optimize the rejection-only control in frequency and pulse width too.
        // This prevents a slightly off-bin local emitter looking more dispersed
        // merely because its sharper zero-DM peak misses the selected window.
        let zero_refined = fold::refine(
            &zero[validation_start..],
            tsamp,
            1. / period,
            resolution * 2.,
            40,
            bins,
        );
        let zero_period = zero_refined.period_s;
        let zero_profile = fold::profile_at(
            &zero[validation_start..],
            tsamp,
            zero_period,
            bins,
            validation_start,
        );
        let zero_window = select_window(&zero_profile);
        let zero_measure = measure(
            &zero[validation_start..],
            validation_start,
            tsamp,
            zero_period,
            &zero_window,
        );
        let contrast = validation_measure.snr.max(0.) / zero_measure.snr.max(1.);
        let dispersion_widths = (*d.iter().max().unwrap() as f64 * tsamp)
            / (period * window.width_bins as f64 / bins as f64);
        let mut time_parts = vec![];
        let mut time_profiles = vec![];
        for part in 0..4 {
            let lo = validation_start + part * validation / 4;
            let hi = validation_start + (part + 1) * validation / 4;
            time_parts.push(measure(&series[lo..hi], lo, tsamp, period, &window));
            time_profiles.push(fold::profile_at(&series[lo..hi], tsamp, period, bins, lo).bins);
        }
        let mut subbands = vec![];
        let mut subband_profiles = vec![];
        for part in 0..4 {
            let lo = part * frequencies.len() / 4;
            let hi = (part + 1) * frequencies.len() / 4;
            let x = dedisperse::apply(
                &data[lo * nsamples..hi * nsamples],
                nsamples,
                &d[lo..hi],
                Some(length),
            )?;
            subbands.push(measure(
                &x[validation_start..],
                validation_start,
                tsamp,
                period,
                &window,
            ));
            subband_profiles.push(
                fold::profile_at(
                    &x[validation_start..],
                    tsamp,
                    period,
                    bins,
                    validation_start,
                )
                .bins,
            );
        }
        // Full-range folded response from first-half per-channel histograms.
        // Channel-delay phase shifts are interpolated; edge truncation is omitted.
        let channel_profiles: Vec<_> = data
            .chunks_exact(nsamples)
            .map(|ch| fold::profile(&ch[..discovery], tsamp, period, bins).bins)
            .collect();
        let response: Vec<_> = trials
            .iter()
            .map(|&dm| {
                let delays = dedisperse::delays(frequencies, tsamp, dm)?;
                let mut combined = vec![0f64; bins];
                for (profile, &delay) in channel_profiles.iter().zip(&delays) {
                    let mean = profile.iter().map(|&v| v as f64).sum::<f64>() / bins as f64;
                    let shift = delay as f64 * tsamp / period * bins as f64;
                    for (b, out) in combined.iter_mut().enumerate() {
                        let pos = (b as f64 + shift).rem_euclid(bins as f64);
                        let index = pos.floor() as usize;
                        let fraction = pos - index as f64;
                        *out += (profile[index] as f64 - mean) * (1. - fraction)
                            + (profile[(index + 1) % bins] as f64 - mean) * fraction;
                    }
                }
                let on = (0..window.width_bins)
                    .map(|i| combined[(window.start_bin + i) % bins])
                    .sum::<f64>()
                    / window.width_bins as f64;
                let off = (combined.iter().sum::<f64>() - on * window.width_bins as f64)
                    / (bins - window.width_bins) as f64;
                let matching = members
                    .iter()
                    .filter(|c| {
                        (c.dm - dm).abs() < step * 0.01
                            && (c.frequency_hz - candidate.frequency_hz).abs() < resolution * 1.5
                    })
                    .max_by(|a, b| a.score.total_cmp(&b.score));
                Ok(DmPoint {
                    dm,
                    spectral_score: matching.map(|c| c.score),
                    spectral_frequency_hz: matching.map(|c| c.frequency_hz),
                    discovery_profile_response: on - off,
                })
            })
            .collect::<Result<_, dedisperse::Error>>()?;
        let mut reasons = vec![];
        if validation_measure.snr < rules.minimum_validation_snr {
            reasons.push("held-out pulse SNR below threshold".into());
        }
        if time_parts
            .iter()
            .filter(|v| v.snr >= rules.minimum_part_snr)
            .count()
            < rules.minimum_time_parts
        {
            reasons.push("pulse does not persist in enough held-out time segments".into());
        }
        if subbands
            .iter()
            .filter(|v| v.snr >= rules.minimum_part_snr)
            .count()
            < rules.minimum_subbands
        {
            reasons.push("pulse does not persist in enough held-out subbands".into());
        }
        if best_dm <= 0. || contrast < rules.minimum_dm_contrast {
            reasons.push("insufficient improvement over zero-DM interference control".into());
        }
        if dispersion_widths < rules.minimum_dispersion_widths {
            reasons.push("dispersion sweep too small relative to selected pulse width to distinguish local interference".into());
        }
        if zero_ratio > rules.maximum_zero_dm_spectral_ratio {
            reasons.push(
                "stronger rationally related periodic signal in zero-DM discovery control".into(),
            );
        }
        results.push(ConfirmedCandidate {
            candidate,
            refined_dm: best_dm,
            profile: best,
            window,
            validation: validation_measure,
            zero_dm_validation: zero_measure,
            zero_dm_period_s: zero_period,
            zero_dm_window: zero_window,
            dm_contrast: contrast,
            dispersion_widths,
            zero_dm_related_spectral_ratio: zero_ratio,
            time_parts,
            subbands,
            time_profiles,
            subband_profiles,
            dm_response: response,
            harmonic_members: members,
            accepted: reasons.is_empty(),
            rejection_reasons: reasons,
        });
    }
    Ok(Outcome {
        rules: rules.clone(),
        discovery_samples: discovery,
        validation_start_sample: validation_start,
        validation_samples: validation,
        common_samples: length,
        trial_count: trials.len(),
        candidates: results,
    })
}

/// Serial reference pipeline, used by evaluation and integration tests. Search
/// input excludes the confirmation half, including maximum channel-delay edges.
pub fn run(
    observation: Observation<'_>,
    trials: &[f64],
    config: &Config,
    limit: usize,
    bins: usize,
    rules: &Rules,
) -> Result<Outcome, dedisperse::Error> {
    let Observation {
        data,
        nsamples,
        frequencies,
        tsamp,
    } = observation;
    let max = trials
        .last()
        .ok_or_else(|| dedisperse::Error("empty DM plan".into()))?;
    let d = dedisperse::delays(frequencies, tsamp, *max)?;
    let length = nsamples
        .checked_sub(*d.iter().max().unwrap())
        .filter(|&n| n >= 512)
        .ok_or_else(|| dedisperse::Error("observation too short".into()))?;
    let mut generated = vec![];
    let discovery = length
        .checked_sub(*d.iter().max().unwrap())
        .filter(|&n| n >= 256)
        .ok_or_else(|| dedisperse::Error("DM delay leaves no independent split".into()))?
        / 2;
    for &dm in trials {
        let x = dedisperse::apply(
            data,
            nsamples,
            &dedisperse::delays(frequencies, tsamp, dm)?,
            Some(discovery),
        )?;
        generated.extend(crate::periodicity(&x, tsamp, dm, config));
    }
    analyze(observation, trials, &generated, limit, bins, rules)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fixed_window_preserves_absolute_phase_and_sign() {
        let window = Window {
            start_bin: 2,
            width_bins: 2,
            bins: 16,
        };
        let x: Vec<_> = (0..10000)
            .map(|i| {
                if window.contains(i % 16) {
                    2.
                } else {
                    (i % 7) as f32 * 0.1
                }
            })
            .collect();
        let a = measure(&x[137..], 137, 0.001, 0.016, &window);
        assert!(a.snr > 50.);
        let b = measure(
            &x[137..],
            137,
            0.001,
            0.016,
            &Window {
                start_bin: 10,
                ..window
            },
        );
        assert!(b.snr < 0.);
    }
    #[test]
    fn grouping_keeps_dm_history_and_harmonics() {
        let mk = |f, dm, score| Candidate {
            frequency_hz: f,
            period_s: 1. / f,
            dm,
            harmonics: 2,
            score,
        };
        let g = group(
            &[
                mk(100., 10., 9.),
                mk(100.01, 12., 10.),
                mk(200., 12., 8.),
                mk(137., 20., 7.),
            ],
            0.1,
            2,
        );
        assert_eq!(g.len(), 2);
        assert_eq!(g[0].1.len(), 3);
    }
}
