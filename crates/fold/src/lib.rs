//! Phase folding and local frequency refinement of dedispersed data.
use serde::Serialize;
#[derive(Clone, Debug, Serialize)]
pub struct Profile {
    pub period_s: f64,
    pub snr: f64,
    /// Approximate reduced chi-square against a flat profile, used for refinement.
    pub reduced_chi_square: f64,
    pub bins: Vec<f32>,
    pub counts: Vec<usize>,
}
pub fn profile(series: &[f32], tsamp: f64, period: f64, bins: usize) -> Profile {
    profile_at(series, tsamp, period, bins, 0)
}
/// Preserve absolute phase when folding a contiguous time segment.
pub fn profile_at(
    series: &[f32],
    tsamp: f64,
    period: f64,
    bins: usize,
    sample_offset: usize,
) -> Profile {
    assert!(bins >= 4 && period.is_finite() && period > 0. && tsamp.is_finite() && tsamp > 0.);
    let mut sums = vec![0.; bins];
    let mut counts = vec![0; bins];
    for (t, &v) in series.iter().enumerate() {
        let phase = ((t + sample_offset) as f64 * tsamp / period).fract();
        let b = ((phase * bins as f64) as usize).min(bins - 1);
        sums[b] += v as f64;
        counts[b] += 1;
    }
    let values: Vec<f32> = sums
        .iter()
        .zip(&counts)
        .map(|(s, &n)| if n > 0 { (s / n as f64) as f32 } else { 0. })
        .collect();
    let center = rfi::median(&values);
    // Lowest 75% of populated bins form an off-pulse baseline. This estimator is
    // for ranking profiles; selection bias means it is not a calibrated significance.
    let mut off: Vec<f32> = values
        .iter()
        .zip(&counts)
        .filter_map(|(&v, &n)| (n > 0).then_some(v))
        .collect();
    off.sort_by(f32::total_cmp);
    off.truncate((off.len() * 3 / 4).max(1));
    let mean = off.iter().map(|&v| v as f64).sum::<f64>() / off.len() as f64;
    let sigma =
        (off.iter().map(|&v| (v as f64 - mean).powi(2)).sum::<f64>() / off.len() as f64).sqrt();
    let peak = values.iter().copied().fold(f32::NEG_INFINITY, f32::max);
    let snr = if sigma > 1e-12 {
        (peak as f64 - center as f64) / sigma
    } else {
        0.
    };
    let sample_mean = series.iter().map(|&v| v as f64).sum::<f64>() / series.len() as f64;
    let sample_variance = series
        .iter()
        .map(|&v| (v as f64 - sample_mean).powi(2))
        .sum::<f64>()
        / series.len() as f64;
    let populated = counts.iter().filter(|&&n| n > 0).count();
    let reduced_chi_square = if sample_variance > 1e-12 && populated > 1 {
        values
            .iter()
            .zip(&counts)
            .map(|(&v, &n)| n as f64 * (v as f64 - sample_mean).powi(2) / sample_variance)
            .sum::<f64>()
            / (populated - 1) as f64
    } else {
        0.0
    };
    Profile {
        period_s: period,
        snr,
        reduced_chi_square,
        bins: values,
        counts,
    }
}
pub fn refine(
    series: &[f32],
    tsamp: f64,
    frequency: f64,
    span_hz: f64,
    steps: usize,
    bins: usize,
) -> Profile {
    assert!(steps > 0 && frequency > span_hz && span_hz >= 0.);
    (0..=steps)
        .map(|i| {
            let f = frequency - span_hz + 2. * span_hz * i as f64 / steps as f64;
            profile(series, tsamp, 1. / f, bins)
        })
        .max_by(|a, b| a.reduced_chi_square.total_cmp(&b.reduced_chi_square))
        .unwrap()
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn true_period_concentrates_pulses() {
        let x: Vec<f32> = (0..20000)
            .map(|i| {
                if i % 100 < 5 {
                    10.
                } else {
                    ((i * 17 % 23) as f32 - 11.) / 10.
                }
            })
            .collect();
        let a = profile(&x, 0.001, 0.1, 32);
        let b = profile(&x, 0.001, 0.1037, 32);
        assert!(a.snr > b.snr * 2.);
        assert_eq!(a.counts.iter().sum::<usize>(), x.len());
    }
}
