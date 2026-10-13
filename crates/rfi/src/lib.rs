//! Conservative robust channel and broadband time-block masks.
use serde::Serialize;

pub fn median(values: &[f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let mut v = values.to_vec();
    let i = v.len() / 2;
    v.select_nth_unstable_by(i, f32::total_cmp);
    v[i]
}
fn outliers(v: &[f32], threshold: f32) -> Vec<bool> {
    let center = median(v);
    let mad = median(&v.iter().map(|x| (x - center).abs()).collect::<Vec<_>>());
    // The relative floor only guards degenerate (e.g. constant or heavily quantized) data. A
    // 1% floor swamped real scatter of block means on high-baseline data and hid broadband bursts.
    let scale = (1.4826 * mad).max(center.abs() * 1e-4).max(1e-6);
    v.iter()
        .map(|x| (x - center).abs() > threshold * scale)
        .collect()
}
#[derive(Debug, Serialize)]
pub struct Mask {
    pub channels: Vec<usize>,
    pub blocks: Vec<usize>,
    pub block_samples: usize,
    pub masked_fraction: f64,
}
/// Replace masked values with each channel's median. Mutates only the supplied array.
pub fn mask(data: &mut [f32], nsamples: usize, block_samples: usize, threshold: f32) -> Mask {
    assert!(nsamples > 0 && data.len().is_multiple_of(nsamples));
    assert!(block_samples > 0 && threshold.is_finite() && threshold > 0.0);
    let nchans = data.len() / nsamples;
    let means: Vec<f32> = data
        .chunks_exact(nsamples)
        .map(|c| c.iter().sum::<f32>() / nsamples as f32)
        .collect();
    let variances: Vec<f32> = data
        .chunks_exact(nsamples)
        .zip(&means)
        .map(|(c, m)| c.iter().map(|v| (v - m).powi(2)).sum::<f32>() / nsamples as f32)
        .collect();
    let means_bad = outliers(&means, threshold);
    let vars_bad = outliers(&variances, threshold);
    let channels: Vec<usize> = (0..nchans)
        .filter(|&c| means_bad[c] || vars_bad[c] || variances[c] < 1e-12)
        .collect();
    let mut bad = vec![false; nchans];
    for &c in &channels {
        bad[c] = true;
    }
    let blocks_n = nsamples.div_ceil(block_samples);
    let active = nchans - channels.len();
    let block_means: Vec<f32> = (0..blocks_n)
        .map(|b| {
            let start = b * block_samples;
            let end = (start + block_samples).min(nsamples);
            let mut sum = 0.0;
            for c in 0..nchans {
                if !bad[c] {
                    sum += data[c * nsamples + start..c * nsamples + end]
                        .iter()
                        .sum::<f32>();
                }
            }
            sum / ((end - start) * active.max(1)) as f32
        })
        .collect();
    let blocks: Vec<usize> = outliers(&block_means, threshold)
        .iter()
        .enumerate()
        .filter_map(|(i, b)| b.then_some(i))
        .collect();
    let mut count = channels.len() * nsamples;
    for c in 0..nchans {
        let row = &mut data[c * nsamples..(c + 1) * nsamples];
        let typical = median(row);
        if bad[c] {
            row.fill(typical);
        } else {
            for &b in &blocks {
                let start = b * block_samples;
                let end = (start + block_samples).min(nsamples);
                row[start..end].fill(typical);
                count += end - start;
            }
        }
    }
    Mask {
        channels,
        blocks,
        block_samples,
        masked_fraction: count as f64 / data.len() as f64,
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn modest_broadband_burst_on_high_baseline_is_masked() {
        // Block mean shifts by +6 on a baseline of 100 with unit noise: below 1% of baseline.
        let (n, chans) = (4000, 16);
        let mut data: Vec<f32> = (0..n * chans)
            .map(|i| 100. + ((i * 2654435761usize) % 1000) as f32 / 500. - 1.)
            .collect();
        for c in 0..chans {
            for t in 1600..1680 {
                data[c * n + t] += 30.;
            }
        }
        let m = mask(&mut data, n, 400, 8.);
        assert!(m.blocks.contains(&4), "{:?}", m.blocks);
    }
    #[test]
    fn planted_channel_and_block_are_masked() {
        let n = 1000;
        let mut data: Vec<f32> = (0..n * 16).map(|i| 10. + (i % 7) as f32).collect();
        data[..n].fill(1000.);
        for c in 1..16 {
            for t in 400..450 {
                data[c * n + t] += 500.;
            }
        }
        let m = mask(&mut data, n, 50, 6.);
        assert!(m.channels.contains(&0));
        assert!(m.blocks.contains(&8));
        assert!(m.masked_fraction > 0.05);
    }
}
