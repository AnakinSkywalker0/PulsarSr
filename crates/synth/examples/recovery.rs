//! Deterministic end-to-end injection evaluation, with noise-only controls.
//! cargo run --release -p synth --example recovery -- reference/recovery.json
use serde::Serialize;
use std::{error::Error, fs::File};

#[derive(Serialize)]
struct Trial {
    seed: u64,
    amplitude: f64,
    injected_period: f64,
    injected_dm: f64,
    duty_cycle: f64,
    expected_matched_filter_snr: f64,
    recovered: bool,
    period_recovered: bool,
    period_error: Option<f64>,
    dm_error: Option<f64>,
    candidates: usize,
    strongest_score: f64,
}
#[derive(Serialize)]
struct Evaluation {
    protocol: String,
    dm_step: f64,
    score_threshold: f64,
    candidate_limit: usize,
    frequency_tolerance_bins: f64,
    dm_tolerance: f64,
    trials: Vec<Trial>,
}
fn main() -> Result<(), Box<dyn Error>> {
    let path = std::env::args()
        .nth(1)
        .unwrap_or("reference/recovery.json".into());
    let mut rows = vec![];
    // Period, DM and duty cycle vary independently across the fixed 40 seeds.
    for amplitude in [0., 0.025, 0.05, 0.1, 0.2, 0.4] {
        for i in 0..40u64 {
            let c = synth::Config {
                nsamples: 16384,
                nchans: 16,
                foff: -8.,
                amplitude,
                seed: 1000 + i,
                period_s: [0.02371, 0.032, 0.04719, 0.06137][i as usize % 4],
                dm: [18., 42., 62., 96., 118.][i as usize % 5],
                duty_cycle: [0.03, 0.05, 0.1][i as usize % 3],
                ..synth::Config::default()
            };
            let mut fb = synth::generate(&c).map_err(std::io::Error::other)?;
            rfi::mask(&mut fb.data, c.nsamples, 400, 8.);
            let frequencies = fb.header.frequencies();
            let max_delay = *dedisperse::delays(&frequencies, c.tsamp, 120.)?
                .iter()
                .max()
                .unwrap();
            let length = c.nsamples - max_delay;
            let resolution = 1. / (length as f64 * c.tsamp);
            let config = search::Config {
                min_hz: 10.,
                max_hz: 1000.,
                threshold: 8.,
                max_candidates: 32,
                birdies: vec![],
            };
            let mut all = vec![];
            for dm in dedisperse::plan(0., 120., 2.)? {
                let d = dedisperse::delays(&frequencies, c.tsamp, dm)?;
                let x = dedisperse::apply(&fb.data, c.nsamples, &d, Some(length))?;
                all.extend(search::periodicity(&x, c.tsamp, dm, &config));
            }
            let candidates = search::sift(all, resolution * 1.5, 10, true);
            let matching = candidates.iter().find(|v| {
                (v.frequency_hz - 1. / c.period_s).abs() <= resolution * 1.5
                    && (v.dm - c.dm).abs() <= 4.
            });
            let period_recovered = candidates
                .iter()
                .any(|v| (v.frequency_hz - 1. / c.period_s).abs() <= resolution * 1.5);
            let period_error = matching.map(|v| v.period_s - c.period_s);
            let dm_error = matching.map(|v| v.dm - c.dm);
            // Gaussian pulse energy in white noise: A sqrt(Nchan Ntime sigma_phase sqrt(pi)).
            // This is an analytic injected-strength proxy, not a measured search significance.
            let expected = amplitude
                * (c.nchans as f64
                    * length as f64
                    * (c.duty_cycle / 2.355)
                    * std::f64::consts::PI.sqrt())
                .sqrt();
            rows.push(Trial {
                seed: c.seed,
                amplitude,
                injected_period: c.period_s,
                injected_dm: c.dm,
                duty_cycle: c.duty_cycle,
                expected_matched_filter_snr: expected,
                recovered: amplitude > 0. && matching.is_some(),
                period_recovered: amplitude > 0. && period_recovered,
                period_error,
                dm_error,
                candidates: candidates.len(),
                strongest_score: candidates.first().map_or(0., |v| v.score),
            });
        }
        let group = &rows[rows.len() - 40..];
        println!(
            "Amplitude {amplitude:.3}: recovered {}/40; noise-only candidate observations {}/40",
            group.iter().filter(|v| v.recovered).count(),
            group
                .iter()
                .filter(|v| v.candidates > 0 && amplitude == 0.)
                .count()
        );
    }
    let evaluation=Evaluation {protocol:"240 paired seeded observations (40 noise-only + 200 injections); fixed settings, development evaluation, Gaussian pulses/noise only; no binary motion or real RFI; strict fundamental + DM matching, no harmonic credit".into(),dm_step:2.,score_threshold:8.,candidate_limit:10,frequency_tolerance_bins:1.5,dm_tolerance:4.,trials:rows};
    serde_json::to_writer_pretty(File::create(path)?, &evaluation)?;
    Ok(())
}
