use clap::{Parser, Subcommand};
use plotters::prelude::*;
use rayon::prelude::*;
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::{
    error::Error,
    fs::{self, File},
    io::{BufWriter, Read, Write},
    path::{Path, PathBuf},
    time::Instant,
};
type Result<T> = std::result::Result<T, Box<dyn Error>>;

#[derive(Parser)]
#[command(
    version,
    about = "Rust filterbank pulsar search and reproducible synthetic fixtures"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}
#[derive(Subcommand)]
enum Command {
    /// Fold a specified period/DM, without searching or declaring a detection.
    Fold {
        file: PathBuf,
        #[arg(long)]
        period: f64,
        #[arg(long)]
        dm: f64,
        #[arg(long, default_value_t = 64)]
        bins: usize,
        #[arg(long, default_value = "output_fold")]
        out: PathBuf,
    },
    /// Print parsed SIGPROC header as JSON.
    Header { file: PathBuf },
    /// Waterfall and DM-zero time-series PNGs.
    Inspect {
        file: PathBuf,
        #[arg(long, default_value = "output")]
        out: PathBuf,
        #[arg(long)]
        max_samples: Option<usize>,
    },
    /// Seeded Gaussian-noise filterbank with one dispersed pulsar.
    Synth {
        file: PathBuf,
        #[arg(long, default_value_t = 65536)]
        samples: usize,
        #[arg(long, default_value_t = 32)]
        channels: usize,
        #[arg(long, default_value_t = 0.00025)]
        tsamp: f64,
        #[arg(long, default_value_t = 0.032)]
        period: f64,
        #[arg(long, default_value_t = 62.)]
        dm: f64,
        #[arg(long, default_value_t = 0.5)]
        amplitude: f64,
        #[arg(long, default_value_t = 0.05)]
        duty: f64,
        #[arg(long, default_value_t = 42)]
        seed: u64,
    },
    /// Dedisperse one DM, exporting little-endian f32 time series and PNG.
    Dedisperse {
        file: PathBuf,
        #[arg(long)]
        dm: f64,
        #[arg(long, default_value = "output")]
        out: PathBuf,
    },
    /// Search trial DMs, sift candidates, refine and fold the strongest candidates.
    Search {
        file: PathBuf,
        #[arg(long, default_value_t = 0.)]
        dm_min: f64,
        #[arg(long, default_value_t = 500.)]
        dm_max: f64,
        /// Omit for sample-time-limited spacing; larger values are explicitly coarser.
        #[arg(long)]
        dm_step: Option<f64>,
        #[arg(long, default_value_t = 15.)]
        min_hz: f64,
        #[arg(long, default_value_t = 2000.)]
        max_hz: f64,
        #[arg(long, default_value_t = 8.)]
        threshold: f64,
        #[arg(long, default_value_t = 5)]
        top: usize,
        #[arg(long, default_value_t = 64)]
        bins: usize,
        #[arg(long)]
        no_rfi: bool,
        /// Spectral masks as CENTER_HZ:HALF_WIDTH_HZ; repeatable.
        #[arg(long,value_parser=parse_birdie)]
        birdie: Vec<(f64, f64)>,
        #[arg(long)]
        max_samples: Option<usize>,
        #[arg(long)]
        threads: Option<usize>,
        #[arg(long, default_value = "output")]
        out: PathBuf,
    },
}
fn parse_birdie(s: &str) -> std::result::Result<(f64, f64), String> {
    let (a, b) = s.split_once(':').ok_or("expected HZ:HALF_WIDTH_HZ")?;
    let a: f64 = a.parse().map_err(|_| "invalid frequency")?;
    let b: f64 = b.parse().map_err(|_| "invalid width")?;
    if !a.is_finite() || !b.is_finite() || a < 0. || b < 0. {
        return Err("nonnegative finite birdie required".into());
    }
    Ok((a, b))
}
fn json(path: impl AsRef<Path>, value: &impl Serialize) -> Result<()> {
    serde_json::to_writer_pretty(BufWriter::new(File::create(path)?), value)?;
    Ok(())
}
fn series(path: impl AsRef<Path>, x: &[f32]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    for v in x {
        w.write_all(&v.to_le_bytes())?;
    }
    w.flush()?;
    Ok(())
}
fn line_plot(
    path: &Path,
    title: &str,
    xlabel: &str,
    ylabel: &str,
    points: &[(f64, f64)],
) -> Result<()> {
    if points.len() < 2 {
        return Ok(());
    }
    let xmin = points.iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
    let xmax = points.iter().map(|p| p.0).fold(f64::NEG_INFINITY, f64::max);
    let ymin = points.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
    let ymax = points.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max);
    let pad = ((ymax - ymin) * 0.1).max(1e-3);
    let root = BitMapBackend::new(path, (1100, 500)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 24))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(65)
        .build_cartesian_2d(xmin..xmax.max(xmin + 1e-6), ymin - pad..ymax + pad)?;
    chart
        .configure_mesh()
        .x_desc(xlabel)
        .y_desc(ylabel)
        .draw()?;
    chart.draw_series(LineSeries::new(points.iter().copied(), &BLUE))?;
    root.present()?;
    Ok(())
}
fn time_plot(path: &Path, x: &[f32], tsamp: f64) -> Result<()> {
    let stride = x.len().div_ceil(4000);
    let points: Vec<_> = x
        .chunks(stride)
        .enumerate()
        .map(|(i, c)| {
            (
                i as f64 * stride as f64 * tsamp,
                c.iter().map(|&v| v as f64).sum::<f64>() / c.len() as f64,
            )
        })
        .collect();
    line_plot(
        path,
        "Dedispersed time series (averaged for display)",
        "Time (s)",
        "Summed power",
        &points,
    )
}
fn inspect(fb: &fbank::Filterbank, out: &Path) -> Result<()> {
    fs::create_dir_all(out)?;
    json(out.join("header.json"), &fb.header)?;
    let dm0 = dedisperse::apply(&fb.data, fb.nsamples, &vec![0; fb.header.nchans], None)?;
    time_plot(&out.join("dm0.png"), &dm0, fb.header.tsamp)?;
    let path = out.join("waterfall.png");
    let root = BitMapBackend::new(&path, (1100, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    let freqs = fb.header.frequencies();
    let low = freqs.iter().copied().fold(f64::INFINITY, f64::min);
    let high = freqs.iter().copied().fold(0., f64::max);
    let duration = fb.nsamples as f64 * fb.header.tsamp;
    let mut chart = ChartBuilder::on(&root)
        .caption("Filterbank waterfall", ("sans-serif", 24))
        .margin(15)
        .x_label_area_size(40)
        .y_label_area_size(65)
        .build_cartesian_2d(0.0..duration, low - 0.5..high + 0.5)?;
    chart
        .configure_mesh()
        .x_desc("Time (s)")
        .y_desc("Frequency (MHz)")
        .draw()?;
    let stride = fb.nsamples.div_ceil(800);
    for (ch, &f) in freqs.iter().enumerate() {
        let row = fb.channel(ch);
        let mean = row.iter().sum::<f32>() / row.len() as f32;
        let sigma = (row.iter().map(|v| (v - mean).powi(2)).sum::<f32>() / row.len() as f32)
            .sqrt()
            .max(1e-6);
        chart.draw_series(row.chunks(stride).enumerate().map(|(i, c)| {
            let z = (c.iter().sum::<f32>() / c.len() as f32 - mean) / sigma;
            let intensity = (0.5 + z as f64 * 0.3).clamp(0., 1.);
            let t0 = i as f64 * stride as f64 * fb.header.tsamp;
            let t1 = (t0 + stride as f64 * fb.header.tsamp).min(duration);
            let half = fb.header.foff.abs().max(1.) / 2.;
            Rectangle::new(
                [(t0, f - half), (t1, f + half)],
                HSLColor(0.65 - intensity * 0.65, 0.9, 0.15 + intensity * 0.55).filled(),
            )
        }))?;
    }
    root.present()?;
    Ok(())
}
fn diagnostic_plot(path: &Path, title: &str, profiles: &[Vec<f32>]) -> Result<()> {
    let root = BitMapBackend::new(path, (1100, 700)).into_drawing_area();
    root.fill(&WHITE)?;
    for (i, (area, profile)) in root
        .split_evenly((profiles.len(), 1))
        .into_iter()
        .zip(profiles)
        .enumerate()
    {
        let mean = profile.iter().map(|&v| v as f64).sum::<f64>() / profile.len() as f64;
        let range = profile
            .iter()
            .map(|&v| (v as f64 - mean).abs())
            .fold(1e-4, f64::max);
        let mut chart = ChartBuilder::on(&area)
            .caption(format!("{title} {}", i + 1), ("sans-serif", 16))
            .margin(8)
            .x_label_area_size(20)
            .y_label_area_size(60)
            .build_cartesian_2d(0.0..1.0, -range..range)?;
        chart.configure_mesh().draw()?;
        chart.draw_series(LineSeries::new(
            profile
                .iter()
                .enumerate()
                .map(|(b, &v)| (b as f64 / profile.len() as f64, v as f64 - mean)),
            &BLUE,
        ))?;
    }
    root.present()?;
    Ok(())
}
#[derive(Serialize)]
struct Folded {
    candidate: search::Candidate,
    refined_dm: f64,
    profile: fold::Profile,
    dm_curve: Vec<(f64, f64)>,
    confirmation: search::confirmation::ConfirmedCandidate,
}
#[derive(Serialize)]
struct Report {
    schema_version: u32,
    input_sha256: String,
    code_revision: Option<String>,
    confirmation: search::confirmation::Outcome,
    input: String,
    version: String,
    sample_count: usize,
    tsamp: f64,
    dm_min: f64,
    dm_max: f64,
    dm_step: f64,
    trial_count: usize,
    common_samples: usize,
    threads: usize,
    elapsed_seconds: f64,
    min_hz: f64,
    max_hz: f64,
    threshold: f64,
    birdies: Vec<(f64, f64)>,
    score_definition: String,
    rfi: Option<rfi::Mask>,
    candidates: Vec<search::Candidate>,
    folded: Vec<Folded>,
}
fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Fold {
            file,
            period,
            dm,
            bins,
            out,
        } => {
            if !period.is_finite() || period <= 0. || bins < 8 {
                return Err("positive period and >=8 bins required".into());
            }
            let fb = fbank::read(&file, None)?;
            let x = dedisperse::apply(
                &fb.data,
                fb.nsamples,
                &dedisperse::delays(&fb.header.frequencies(), fb.header.tsamp, dm)?,
                None,
            )?;
            let p = fold::profile(&x, fb.header.tsamp, period, bins);
            fs::create_dir_all(&out)?;
            json(
                out.join("fold.json"),
                &serde_json::json!({"input":file,"dm":dm,"profile":p,"searched":false}),
            )?;
            let points: Vec<_> = p
                .bins
                .iter()
                .enumerate()
                .map(|(i, &v)| (i as f64 / bins as f64, v as f64))
                .collect();
            line_plot(
                &out.join("profile.png"),
                "Specified period/DM fold",
                "Pulse phase",
                "Mean power",
                &points,
            )?;
        }
        Command::Header { file } => println!(
            "{}",
            serde_json::to_string_pretty(&fbank::read_header(file)?)?
        ),
        Command::Inspect {
            file,
            out,
            max_samples,
        } => inspect(&fbank::read(file, max_samples)?, &out)?,
        Command::Synth {
            file,
            samples,
            channels,
            tsamp,
            period,
            dm,
            amplitude,
            duty,
            seed,
        } => {
            let config = synth::Config {
                nsamples: samples,
                nchans: channels,
                tsamp,
                period_s: period,
                dm,
                amplitude,
                duty_cycle: duty,
                seed,
                ..synth::Config::default()
            };
            let fb = synth::generate(&config).map_err(std::io::Error::other)?;
            fbank::write(BufWriter::new(File::create(&file)?), &fb)?;
            json(file.with_extension("truth.json"), &config)?;
            println!("Wrote {} and truth metadata", file.display());
        }
        Command::Dedisperse { file, dm, out } => {
            fs::create_dir_all(&out)?;
            let fb = fbank::read(file, None)?;
            let delays = dedisperse::delays(&fb.header.frequencies(), fb.header.tsamp, dm)?;
            let x = dedisperse::apply(&fb.data, fb.nsamples, &delays, None)?;
            series(out.join("series.f32"), &x)?;
            time_plot(&out.join("series.png"), &x, fb.header.tsamp)?;
            json(
                out.join("series.json"),
                &serde_json::json!({"dm":dm,"tsamp":fb.header.tsamp,"samples":x.len(),"dtype":"little-endian float32","normalization":"raw channel sum","delays":delays}),
            )?;
        }
        Command::Search {
            file,
            dm_min,
            dm_max,
            dm_step,
            min_hz,
            max_hz,
            threshold,
            top,
            bins,
            no_rfi,
            birdie,
            max_samples,
            threads,
            out,
        } => {
            if !min_hz.is_finite()
                || !max_hz.is_finite()
                || min_hz <= 0.
                || max_hz <= min_hz
                || !threshold.is_finite()
                || threshold <= 0.
                || top == 0
                || bins < 8
            {
                return Err("invalid search settings".into());
            }
            if let Some(n) = threads {
                if n == 0 {
                    return Err("threads must be positive".into());
                }
                rayon::ThreadPoolBuilder::new()
                    .num_threads(n)
                    .build_global()?;
            }
            let start = Instant::now();
            fs::create_dir_all(&out)?;
            let mut fb = fbank::read(&file, max_samples)?;
            let frequencies = fb.header.frequencies();
            let tsamp = fb.header.tsamp;
            json(out.join("header.json"), &fb.header)?;
            let mask = if no_rfi {
                None
            } else {
                Some(rfi::mask(
                    &mut fb.data,
                    fb.nsamples,
                    ((0.1 / tsamp) as usize).max(32),
                    8.,
                ))
            };
            let step = dm_step.unwrap_or(dedisperse::dm_step(&frequencies, tsamp)?);
            let trials = dedisperse::plan(dm_min, dm_max, step)?;
            let max_delay = *dedisperse::delays(&frequencies, tsamp, dm_max)?
                .iter()
                .max()
                .unwrap();
            let length = fb
                .nsamples
                .checked_sub(max_delay)
                .filter(|&n| n >= 512)
                .ok_or("DM range leaves fewer than 512 samples")?;
            let discovery = length
                .checked_sub(max_delay)
                .filter(|&n| n >= 256)
                .ok_or("DM range leaves no independent discovery/confirmation split")?
                / 2;
            let config = search::Config {
                min_hz,
                max_hz,
                threshold,
                max_candidates: 32,
                birdies: birdie.clone(),
            };
            println!(
                "{} samples, {} channels, {} DM trials (step {:.6}), {} threads",
                fb.nsamples,
                fb.header.nchans,
                trials.len(),
                step,
                rayon::current_num_threads()
            );
            let lists: Vec<_> = trials
                .par_iter()
                .map(|&dm| {
                    let delay = dedisperse::delays(&frequencies, tsamp, dm)?;
                    let x = dedisperse::apply(&fb.data, fb.nsamples, &delay, Some(discovery))?;
                    Ok::<_, dedisperse::Error>(search::periodicity(&x, tsamp, dm, &config))
                })
                .collect();
            let mut all = vec![];
            for list in lists {
                all.extend(list?);
            }
            let confirmation = search::confirmation::analyze(
                search::confirmation::Observation {
                    data: &fb.data,
                    nsamples: fb.nsamples,
                    frequencies: &frequencies,
                    tsamp,
                },
                &trials,
                &all,
                top,
                bins,
                &search::confirmation::Rules::default(),
            )?;
            let candidates: Vec<_> = confirmation
                .candidates
                .iter()
                .map(|c| c.candidate.clone())
                .collect();
            let mut folded = vec![];
            for (index, confirmed) in confirmation.candidates.iter().enumerate() {
                let candidate = &confirmed.candidate;
                let best = &confirmed.profile;
                let best_dm = confirmed.refined_dm;
                let curve: Vec<_> = confirmed
                    .dm_response
                    .iter()
                    .map(|p| (p.dm, p.discovery_profile_response))
                    .collect();
                let phase: Vec<_> = best
                    .bins
                    .iter()
                    .enumerate()
                    .map(|(i, &v)| (i as f64 / bins as f64, v as f64))
                    .collect();
                line_plot(
                    &out.join(format!("candidate_{index}_profile.png")),
                    &format!("P = {:.6} ms, DM = {:.3}", best.period_s * 1000., best_dm),
                    "Pulse phase",
                    "Mean power",
                    &phase,
                )?;
                line_plot(
                    &out.join(format!("candidate_{index}_dm.png")),
                    "Full-range discovery DM response (phase-histogram approximation)",
                    "DM (pc cm^-3)",
                    "Fixed-window power contrast",
                    &curve,
                )?;
                println!(
                    "#{index}: P={:.6} ms DM={:.3} spectral score={:.2} held-out SNR={:.2} {}",
                    best.period_s * 1000.,
                    best_dm,
                    candidate.score,
                    confirmed.validation.snr,
                    if confirmed.accepted {
                        "ACCEPTED for follow-up"
                    } else {
                        "REJECTED"
                    }
                );
                diagnostic_plot(
                    &out.join(format!("candidate_{index}_time.png")),
                    "Held-out time segments",
                    &confirmed.time_profiles,
                )?;
                diagnostic_plot(
                    &out.join(format!("candidate_{index}_subbands.png")),
                    "Held-out frequency subbands",
                    &confirmed.subband_profiles,
                )?;
                folded.push(Folded {
                    candidate: candidate.clone(),
                    refined_dm: best_dm,
                    profile: best.clone(),
                    dm_curve: curve,
                    confirmation: confirmed.clone(),
                });
            }
            let mut reader = File::open(&file)?;
            let mut hash = Sha256::new();
            let mut buffer = [0; 65536];
            loop {
                let n = reader.read(&mut buffer)?;
                if n == 0 {
                    break;
                }
                hash.update(&buffer[..n]);
            }
            let code_revision = std::process::Command::new("git")
                .args(["rev-parse", "HEAD"])
                .output()
                .ok()
                .filter(|o| o.status.success())
                .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string());
            let report=Report {schema_version:2,input_sha256:format!("{:x}",hash.finalize()),code_revision,confirmation,input:file.display().to_string(),version:env!("CARGO_PKG_VERSION").into(),sample_count:fb.nsamples,tsamp,dm_min,dm_max,dm_step:step,trial_count:trials.len(),common_samples:length,threads:rayon::current_num_threads(),elapsed_seconds:start.elapsed().as_secs_f64(),min_hz,max_hz,threshold,birdies:birdie,score_definition:"discovery spectral ranking and approximate profile statistics; confirmation fixed-window on/off mean difference / sample standard error on disjoint second half; colored noise/RFI may violate independent-noise assumptions".into(),rfi:mask,candidates,folded};
            json(out.join("report.json"), &report)?;
            println!("Saved {} in {:.2}s", out.display(), report.elapsed_seconds);
        }
    }
    Ok(())
}
