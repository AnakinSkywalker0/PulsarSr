//! Incoherent integer-sample dedispersion; frequencies are MHz, times seconds.
use thiserror::Error;
pub const DISPERSION_SECONDS_MHZ: f64 = 4148.808;

#[derive(Debug, Error)]
#[error("invalid dedispersion input: {0}")]
pub struct Error(pub String);

pub fn delays(frequencies: &[f64], tsamp: f64, dm: f64) -> Result<Vec<usize>, Error> {
    if frequencies.is_empty()
        || frequencies.iter().any(|f| !f.is_finite() || *f <= 0.0)
        || !tsamp.is_finite()
        || tsamp <= 0.0
        || !dm.is_finite()
        || dm < 0.0
    {
        return Err(Error(
            "positive frequencies/tsamp and nonnegative DM required".into(),
        ));
    }
    let reference = frequencies.iter().copied().fold(0.0, f64::max);
    frequencies
        .iter()
        .map(|f| {
            let samples =
                (DISPERSION_SECONDS_MHZ * dm * (f.powi(-2) - reference.powi(-2)) / tsamp).round();
            if samples > usize::MAX as f64 {
                Err(Error("delay overflow".into()))
            } else {
                Ok(samples.max(0.0) as usize)
            }
        })
        .collect()
}
pub fn dm_step(frequencies: &[f64], tsamp: f64) -> Result<f64, Error> {
    delays(frequencies, tsamp, 0.0)?;
    let lo = frequencies.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = frequencies.iter().copied().fold(0.0, f64::max);
    let span = DISPERSION_SECONDS_MHZ * (lo.powi(-2) - hi.powi(-2));
    if span <= 0.0 {
        return Err(Error("DM planning requires finite bandwidth".into()));
    }
    Ok(tsamp / span)
}
pub fn plan(min: f64, max: f64, step: f64) -> Result<Vec<f64>, Error> {
    if !min.is_finite()
        || !max.is_finite()
        || !step.is_finite()
        || min < 0.0
        || max < min
        || step <= 0.0
    {
        return Err(Error("invalid DM bounds/step".into()));
    }
    let n = ((max - min) / step).floor() as usize + 1;
    if n > 100_000 {
        return Err(Error("more than 100000 DM trials".into()));
    }
    Ok((0..n).map(|i| min + i as f64 * step).collect())
}
/// Common output length can be selected across trials to keep FFT bins identical.
pub fn apply(
    data: &[f32],
    nsamples: usize,
    delays: &[usize],
    length: Option<usize>,
) -> Result<Vec<f32>, Error> {
    if delays.is_empty() || nsamples == 0 || nsamples.checked_mul(delays.len()) != Some(data.len())
    {
        return Err(Error("array shape mismatch".into()));
    }
    let max_delay = *delays.iter().max().unwrap();
    let available = nsamples
        .checked_sub(max_delay)
        .ok_or_else(|| Error("delay exceeds observation".into()))?;
    let n = length.unwrap_or(available);
    if n < 2 || n > available {
        return Err(Error("invalid output length".into()));
    }
    let mut output = vec![0.0; n];
    for (ch, &delay) in delays.iter().enumerate() {
        for (out, &v) in output
            .iter_mut()
            .zip(&data[ch * nsamples + delay..ch * nsamples + delay + n])
        {
            *out += v;
        }
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn dispersed_delta_collapses_for_both_frequency_orders() {
        for freqs in [vec![1400., 1300., 1200.], vec![1200., 1300., 1400.]] {
            let d = delays(&freqs, 0.001, 20.).unwrap();
            let n = 100;
            let mut data = vec![0.; n * 3];
            for (c, &delay) in d.iter().enumerate() {
                data[c * n + 10 + delay] = 1.;
            }
            let out = apply(&data, n, &d, None).unwrap();
            assert_eq!(out[10], 3.);
            assert_eq!(out.iter().sum::<f32>(), 3.);
        }
    }
    #[test]
    fn units_and_invalid_inputs() {
        assert_eq!(delays(&[1000., 2000.], 0.001, 1.).unwrap(), vec![3, 0]);
        assert!(plan(0., 5., 0.).is_err());
        assert!(delays(&[0.], 1., 1.).is_err());
        assert!(apply(&[1., 2.], 2, &[3], None).is_err());
    }
}
