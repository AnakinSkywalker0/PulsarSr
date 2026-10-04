//! Seeded Gaussian noise plus dispersed periodic Gaussian pulses.
use fbank::{Filterbank, Header, Value};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rand_distr::StandardNormal;
use serde::Serialize;
use std::collections::BTreeMap;

#[derive(Clone, Debug, Serialize)]
pub struct Config {
    pub nsamples: usize,
    pub nchans: usize,
    pub tsamp: f64,
    pub fch1: f64,
    pub foff: f64,
    pub period_s: f64,
    pub dm: f64,
    pub duty_cycle: f64,
    pub amplitude: f64,
    pub seed: u64,
}
impl Default for Config {
    fn default() -> Self {
        Self {
            nsamples: 65536,
            nchans: 32,
            tsamp: 0.00025,
            fch1: 1450.,
            foff: -4.,
            period_s: 0.032,
            dm: 62.,
            duty_cycle: 0.05,
            amplitude: 0.5,
            seed: 42,
        }
    }
}
pub fn generate(c: &Config) -> Result<Filterbank, String> {
    if c.nsamples < 2
        || c.nchans < 1
        || !c.period_s.is_finite()
        || c.period_s <= 0.
        || !c.duty_cycle.is_finite()
        || c.duty_cycle <= 0.
        || c.duty_cycle >= 0.5
        || !c.amplitude.is_finite()
        || c.amplitude < 0.
    {
        return Err("invalid synthetic configuration".into());
    }
    let frequencies: Vec<f64> = (0..c.nchans).map(|i| c.fch1 + i as f64 * c.foff).collect();
    dedisperse::delays(&frequencies, c.tsamp, c.dm).map_err(|e| e.to_string())?;
    let reference = frequencies.iter().copied().fold(0., f64::max);
    let mut rng = ChaCha8Rng::seed_from_u64(c.seed);
    let mut data = Vec::with_capacity(c.nsamples * c.nchans);
    for f in frequencies {
        // Fractional physical delays deliberately differ from integer dedispersion.
        let delay = dedisperse::DISPERSION_SECONDS_MHZ * c.dm * (f.powi(-2) - reference.powi(-2));
        for t in 0..c.nsamples {
            let phase = ((t as f64 * c.tsamp - delay) / c.period_s - 0.25).rem_euclid(1.);
            let distance = phase.min(1. - phase);
            let pulse = c.amplitude * (-0.5 * (distance / (c.duty_cycle / 2.355)).powi(2)).exp();
            let noise: f64 = rng.sample(StandardNormal);
            data.push((100. + noise + pulse) as f32);
        }
    }
    let mut fields = BTreeMap::new();
    fields.insert(
        "source_name".into(),
        Value::Text("seeded synthetic pulsar".into()),
    );
    Ok(Filterbank {
        header: Header {
            fields,
            nchans: c.nchans,
            nbits: 32,
            tsamp: c.tsamp,
            fch1: c.fch1,
            foff: c.foff,
            data_offset: 0,
        },
        nsamples: c.nsamples,
        data,
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn full_pipeline_recovers_seeded_injection() {
        let c = Config {
            nsamples: 16384,
            nchans: 16,
            amplitude: 2.,
            ..Config::default()
        };
        let fb = generate(&c).unwrap();
        let d = dedisperse::delays(&fb.header.frequencies(), c.tsamp, c.dm).unwrap();
        let x = dedisperse::apply(&fb.data, c.nsamples, &d, None).unwrap();
        let candidates = search::periodicity(&x, c.tsamp, c.dm, &search::Config::default());
        assert!(
            candidates
                .iter()
                .any(|v| (v.period_s - c.period_s).abs() < 0.0001)
        );
        assert!(fold::profile(&x, c.tsamp, c.period_s, 32).snr > 10.);
        let mut bytes = vec![];
        fbank::write(&mut bytes, &fb).unwrap();
        assert_eq!(fbank::decode(&bytes, None).unwrap().data, fb.data);
    }
    #[test]
    fn confirmation_accepts_pulsar_rejects_noise_and_zero_dm_emitter() {
        let c = Config {
            nsamples: 16384,
            nchans: 16,
            foff: -8.,
            amplitude: 0.8,
            seed: 7812,
            ..Config::default()
        };
        let trials = dedisperse::plan(0., 120., 2.).unwrap();
        let config = search::Config {
            min_hz: 10.,
            max_hz: 1000.,
            threshold: 8.,
            ..search::Config::default()
        };
        for (dm, amplitude, expect_signal) in [(62., 0.8, true), (0., 0.8, false), (62., 0., false)]
        {
            let fb = generate(&Config {
                dm,
                amplitude,
                ..c.clone()
            })
            .unwrap();
            let frequencies = fb.header.frequencies();
            let outcome = search::confirmation::run(
                search::confirmation::Observation {
                    data: &fb.data,
                    nsamples: fb.nsamples,
                    frequencies: &frequencies,
                    tsamp: c.tsamp,
                },
                &trials,
                &config,
                10,
                32,
                &search::confirmation::Rules::default(),
            )
            .unwrap();
            assert_eq!(
                outcome.candidates.iter().any(|v| v.accepted),
                expect_signal,
                "DM {dm} amplitude {amplitude}"
            );
            assert!(outcome.validation_start_sample > outcome.discovery_samples);
            assert!(
                outcome
                    .candidates
                    .iter()
                    .all(|v| v.dm_response.len() == trials.len())
            );
        }
    }
}
