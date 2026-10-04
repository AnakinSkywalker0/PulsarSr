//! FFT periodicity search. Scores are ranking statistics, NOT calibrated Gaussian sigmas.
use realfft::RealFftPlanner;
use serde::Serialize;
pub mod confirmation;

#[derive(Clone, Debug, Serialize)]
pub struct Candidate {
    pub frequency_hz: f64,
    pub period_s: f64,
    pub dm: f64,
    pub harmonics: usize,
    pub score: f64,
}
#[derive(Clone, Debug)]
pub struct Config {
    pub min_hz: f64,
    pub max_hz: f64,
    pub threshold: f64,
    pub max_candidates: usize,
    pub birdies: Vec<(f64, f64)>,
}
impl Default for Config {
    fn default() -> Self {
        Self {
            min_hz: 1.,
            max_hz: f64::INFINITY,
            threshold: 6.,
            max_candidates: 32,
            birdies: vec![],
        }
    }
}
pub fn periodicity(series: &[f32], tsamp: f64, dm: f64, config: &Config) -> Vec<Candidate> {
    assert!(series.len() >= 2 && tsamp.is_finite() && tsamp > 0.0);
    let size = series.len().next_power_of_two() * 2;
    let mean = series.iter().map(|&v| v as f64).sum::<f64>() / series.len() as f64;
    let mut input = vec![0.; size];
    for (dst, &src) in input.iter_mut().zip(series) {
        *dst = src - mean as f32;
    }
    let fft = RealFftPlanner::<f32>::new().plan_fft_forward(size);
    let mut output = fft.make_output_vec();
    fft.process(&mut input, &mut output).expect("FFT shape");
    let mut power: Vec<f32> = output.iter().map(|c| c.norm_sqr()).collect();
    // Sliding-window median estimated every 64 bins; suppress slowly varying red noise.
    let raw = power.clone();
    // Degenerate noise-free fixtures have zero medians. A global floor prevents
    // tiny FFT roundoff/leakage from becoming arbitrarily large normalized peaks.
    let floor = (raw.iter().map(|&p| p as f64).sum::<f64>() / raw.len() as f64 * 1e-6) as f32;
    for (block, slice) in power.chunks_mut(64).enumerate() {
        let center = block * 64;
        let lo = center.saturating_sub(128).max(1);
        let hi = (center + 192).min(raw.len());
        let scale = (rfi::median(&raw[lo..hi]) / std::f32::consts::LN_2)
            .max(floor)
            .max(1e-20);
        for p in slice {
            *p /= scale;
        }
    }
    power[0] = 0.;
    let bin_hz = 1. / (size as f64 * tsamp);
    for (i, p) in power.iter_mut().enumerate() {
        let f = i as f64 * bin_hz;
        if config
            .birdies
            .iter()
            .any(|&(center, width)| (f - center).abs() <= width)
        {
            *p = 0.;
        }
    }
    let first = (config.min_hz / bin_hz).ceil().max(1.) as usize;
    let last = ((config.max_hz.min(0.5 / tsamp) / bin_hz).floor() as usize).min(power.len() - 2);
    let mut candidates = vec![];
    for h in [1, 2, 4, 8, 16] {
        let end = last.min((power.len() - 2) / h);
        if end <= first + 1 {
            continue;
        }
        let scores: Vec<f64> = (first..=end)
            .map(|k| {
                // Require evidence at the proposed fundamental. Otherwise one
                // bright harmonic creates spurious fractional-frequency families.
                if h > 1 && power[k] < 4.0 {
                    return f64::NEG_INFINITY;
                }
                // Zero padding oversamples by >=2. Interpolate the fundamental
                // peak before sampling its harmonics; k*j alone drifts off-bin
                // by j times the fundamental rounding error.
                let curvature = raw[k - 1] as f64 - 2. * raw[k] as f64 + raw[k + 1] as f64;
                let offset = if curvature < 0. {
                    (0.5 * (raw[k - 1] as f64 - raw[k + 1] as f64) / curvature).clamp(-0.5, 0.5)
                } else {
                    0.
                };
                let sum: f64 = (1..=h)
                    .map(|j| {
                        let position = (k as f64 + offset) * j as f64;
                        let index = position.round() as usize;
                        power[index.min(power.len() - 1)] as f64
                    })
                    .sum();
                (sum - h as f64) / (h as f64).sqrt()
            })
            .collect();
        for i in 1..scores.len() - 1 {
            if scores[i] >= config.threshold
                && scores[i] > scores[i - 1]
                && scores[i] >= scores[i + 1]
            {
                let k = first + i;
                let curvature = raw[k - 1] as f64 - 2. * raw[k] as f64 + raw[k + 1] as f64;
                let offset = if curvature < 0. {
                    (0.5 * (raw[k - 1] as f64 - raw[k + 1] as f64) / curvature).clamp(-0.5, 0.5)
                } else {
                    0.
                };
                let frequency = (k as f64 + offset) * bin_hz;
                candidates.push(Candidate {
                    frequency_hz: frequency,
                    period_s: 1. / frequency,
                    dm,
                    harmonics: h,
                    score: scores[i],
                });
            }
        }
    }
    candidates.sort_by(|a, b| b.score.total_cmp(&a.score));
    let mut kept: Vec<Candidate> = vec![];
    for c in candidates {
        if kept
            .iter()
            .any(|k| (k.frequency_hz - c.frequency_hz).abs() < 1. / (series.len() as f64 * tsamp))
        {
            continue;
        }
        kept.push(c);
        if kept.len() >= config.max_candidates {
            break;
        }
    }
    kept
}
/// Keep best representative per frequency/harmonic family across all DM trials.
pub fn sift(
    mut candidates: Vec<Candidate>,
    resolution_hz: f64,
    limit: usize,
    reject_dm_zero: bool,
) -> Vec<Candidate> {
    candidates.sort_by(|a, b| b.score.total_cmp(&a.score));
    let mut kept: Vec<Candidate> = vec![];
    let mut rejected: Vec<f64> = vec![];
    for c in candidates {
        let related = |f: f64| {
            (1..=16).any(|h| {
                (c.frequency_hz - f * h as f64).abs() <= resolution_hz
                    || (f - c.frequency_hz * h as f64).abs() <= resolution_hz
            })
        };
        if kept.iter().any(|k| related(k.frequency_hz)) || rejected.iter().any(|&f| related(f)) {
            continue;
        }
        if reject_dm_zero && c.dm <= 0.0 {
            rejected.push(c.frequency_hz);
            continue;
        }
        kept.push(c);
        if kept.len() >= limit {
            break;
        }
    }
    kept
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn finds_sine_and_respects_birdies() {
        let tsamp = 0.001;
        let f = 31.25;
        let x: Vec<f32> = (0..8192)
            .map(|i| (std::f64::consts::TAU * f * i as f64 * tsamp).sin() as f32)
            .collect();
        let config = Config {
            min_hz: 10.,
            ..Config::default()
        };
        let c = periodicity(&x, tsamp, 0., &config);
        assert!((c[0].frequency_hz - f).abs() < 0.1, "{c:?}");
        let c = periodicity(
            &x,
            tsamp,
            0.,
            &Config {
                birdies: vec![(f, 1.)],
                ..config
            },
        );
        assert!(c.iter().all(|c| (c.frequency_hz - f).abs() > 1.));
    }
    #[test]
    fn sifts_dm_duplicates_and_harmonics() {
        let mk = |f, dm, score| Candidate {
            frequency_hz: f,
            period_s: 1. / f,
            dm,
            harmonics: 1,
            score,
        };
        let c = sift(
            vec![
                mk(100., 20., 10.),
                mk(200., 21., 8.),
                mk(100.01, 19., 9.),
                mk(137., 0., 20.),
            ],
            0.1,
            10,
            true,
        );
        assert_eq!(c.len(), 1);
        assert_eq!(c[0].dm, 20.);
    }

    #[test]
    fn off_bin_pulse_train_retains_fundamental() {
        let dt = 0.00025;
        let period = 0.02371;
        let mut state = 123456789u64;
        let series: Vec<f32> = (0..16384)
            .map(|i| {
                let phase = (i as f64 * dt / period).fract();
                let distance = (phase - 0.25).abs();
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                let noise = (state as f64 / u64::MAX as f64 - 0.5) as f32;
                (-0.5 * (distance / 0.013).powi(2)).exp() as f32 + noise * 0.4
            })
            .collect();
        let candidates = periodicity(&series, dt, 18., &Config::default());
        let kept = sift(candidates, 1.5 / (series.len() as f64 * dt), 10, false);
        assert!((kept[0].frequency_hz - 1. / period).abs() < 0.2, "{kept:?}");
    }
}
