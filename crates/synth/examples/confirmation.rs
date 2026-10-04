//! Fixed split-confirm-v1 rules evaluated once per named split. Development seeds
//! 20000..; held-out seeds 40000.. never used for threshold selection.
use serde::Serialize;
use std::{error::Error, fs::File};

#[derive(Serialize)]
struct Row {
    seed: u64,
    scenario: String,
    amplitude: f64,
    period_s: f64,
    dm: f64,
    duty_cycle: f64,
    expected_matched_filter_snr: f64,
    accepted_candidates: usize,
    recovered: bool,
    period_error_s: Option<f64>,
    dm_error: Option<f64>,
    false_accepted: usize,
    diagnostics: Vec<search::confirmation::ConfirmedCandidate>,
}
#[derive(Serialize)]
struct Evaluation {
    split: String,
    protocol: String,
    rules: search::confirmation::Rules,
    generated_candidates_per_dm: usize,
    confirmation_budget: usize,
    dm_step: f64,
    bins: usize,
    observations: Vec<Row>,
}
fn contaminate(fb: &mut fbank::Filterbank, scenario: &str) {
    let n = fb.nsamples;
    let channels = fb.header.nchans;
    let dt = fb.header.tsamp;
    match scenario {
        "colored-noise" | "pulsar-colored" => {
            for row in fb.data.chunks_exact_mut(n) {
                let mut state = 0.;
                for v in row {
                    let white = *v - 100.;
                    state = 0.95 * state + 0.31225 * white;
                    *v += state;
                }
            }
        }
        "broadband-rfi" => {
            for row in fb.data.chunks_exact_mut(n) {
                for (t, v) in row.iter_mut().enumerate() {
                    let phase = (t as f64 * dt / 0.02731).fract();
                    *v += if phase < 0.05 { 4. } else { 0. };
                }
            }
        }
        "narrowband-rfi" => {
            for t in 0..n {
                fb.data[2 * n + t] +=
                    8. * (std::f64::consts::TAU * 37.21 * t as f64 * dt).sin() as f32;
            }
        }
        "burst-rfi" | "pulsar-rfi" => {
            for c in 0..channels {
                for t in n / 3..n / 3 + 80 {
                    fb.data[c * n + t] += 30.;
                }
            }
            for t in 0..n {
                fb.data[2 * n + t] +=
                    10. * (std::f64::consts::TAU * 51.73 * t as f64 * dt).sin() as f32;
            }
        }
        _ => {}
    }
}
fn main() -> Result<(), Box<dyn Error>> {
    let split = std::env::args().nth(1).unwrap_or("development".into());
    let output = std::env::args()
        .nth(2)
        .unwrap_or(format!("reference/confirmation_{split}.json"));
    let base = match split.as_str() {
        "development" => 20000,
        "held-out" => 40000,
        _ => return Err("split must be development or held-out".into()),
    };
    let rules = search::confirmation::Rules::default();
    let config = search::Config {
        min_hz: 10.,
        max_hz: 1000.,
        threshold: 8.,
        max_candidates: 32,
        birdies: vec![],
    };
    let trials = dedisperse::plan(0., 120., 2.)?;
    let mut rows = vec![];
    for (s, scenario) in [
        "white-noise",
        "colored-noise",
        "broadband-rfi",
        "narrowband-rfi",
        "burst-rfi",
        "pulsar",
        "pulsar-colored",
        "pulsar-rfi",
    ]
    .iter()
    .enumerate()
    {
        for i in 0..32 {
            let positive = scenario.starts_with("pulsar");
            let amplitude = if positive {
                [0.1, 0.2, 0.4, 0.8][i / 8]
            } else {
                0.
            };
            let c = synth::Config {
                nsamples: 16384,
                nchans: 16,
                foff: -8.,
                amplitude,
                seed: base + s as u64 * 100 + i as u64,
                period_s: [0.02371, 0.032, 0.04719, 0.06137][i % 4],
                dm: [18., 42., 62., 96.][i % 4],
                duty_cycle: [0.03, 0.05, 0.1][i % 3],
                ..synth::Config::default()
            };
            let mut fb = synth::generate(&c).map_err(std::io::Error::other)?;
            contaminate(&mut fb, scenario);
            rfi::mask(&mut fb.data, c.nsamples, 400, 8.);
            let result = search::confirmation::run(
                search::confirmation::Observation {
                    data: &fb.data,
                    nsamples: fb.nsamples,
                    frequencies: &fb.header.frequencies(),
                    tsamp: c.tsamp,
                },
                &trials,
                &config,
                10,
                32,
                &rules,
            )?;
            let matched = result.candidates.iter().filter(|v| v.accepted).find(|v| {
                (1. / v.profile.period_s - 1. / c.period_s).abs()
                    <= 1.5 / (result.discovery_samples as f64 * c.tsamp)
                    && (v.refined_dm - c.dm).abs() <= 6.
            });
            let accepted = result.candidates.iter().filter(|v| v.accepted).count();
            let recovered = positive && matched.is_some();
            rows.push(Row {
                seed: c.seed,
                scenario: (*scenario).into(),
                amplitude,
                period_s: c.period_s,
                dm: c.dm,
                duty_cycle: c.duty_cycle,
                expected_matched_filter_snr: amplitude
                    * (c.nchans as f64 * result.validation_samples as f64 * c.duty_cycle / 2.355
                        * std::f64::consts::PI.sqrt())
                    .sqrt(),
                accepted_candidates: accepted,
                recovered,
                period_error_s: matched.map(|v| v.profile.period_s - c.period_s),
                dm_error: matched.map(|v| v.refined_dm - c.dm),
                false_accepted: accepted - usize::from(recovered),
                diagnostics: result.candidates,
            });
        }
        let group = &rows[rows.len() - 32..];
        println!(
            "{split} {scenario}: recovery {}/32; observations with false accepted {}/32",
            group.iter().filter(|r| r.recovered).count(),
            group.iter().filter(|r| r.false_accepted > 0).count()
        );
    }
    let evaluation=Evaluation{split,protocol:"fixed split-confirm-v1; first half discovers/refines/selects fixed window; max-channel-delay guard gap; second half tests, 4 time parts/4 subbands; full-observation unsupervised RFI masks; false alarm unit observation; no rule changes after held-out evaluation".into(),rules,generated_candidates_per_dm:32,confirmation_budget:10,dm_step:2.,bins:32,observations:rows};
    serde_json::to_writer_pretty(File::create(output)?, &evaluation)?;
    Ok(())
}
