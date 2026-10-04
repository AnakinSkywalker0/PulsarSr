//! Strict SIGPROC reader and channel-major arrays for numerical stages.
use memmap2::Mmap;
use serde::{Deserialize, Serialize};
use std::{
    collections::BTreeMap,
    fs::File,
    io::{Cursor, Read, Write},
    path::Path,
};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("{0}")]
    Io(#[from] std::io::Error),
    #[error("invalid or unsupported filterbank: {0}")]
    Format(String),
}
pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(untagged)]
pub enum Value {
    Int(i32),
    Float(f64),
    Text(String),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Header {
    pub fields: BTreeMap<String, Value>,
    pub nchans: usize,
    pub nbits: u32,
    pub tsamp: f64,
    pub fch1: f64,
    pub foff: f64,
    pub data_offset: usize,
}
impl Header {
    pub fn frequencies(&self) -> Vec<f64> {
        (0..self.nchans)
            .map(|c| self.fch1 + c as f64 * self.foff)
            .collect()
    }
}

/// Each contiguous row is one channel: `data[channel * nsamples + time]`.
#[derive(Debug, Clone)]
pub struct Filterbank {
    pub header: Header,
    pub nsamples: usize,
    pub data: Vec<f32>,
}
impl Filterbank {
    pub fn channel(&self, c: usize) -> &[f32] {
        &self.data[c * self.nsamples..(c + 1) * self.nsamples]
    }
}

fn string(c: &mut Cursor<&[u8]>) -> Result<String> {
    let mut b = [0; 4];
    c.read_exact(&mut b)?;
    let n = i32::from_le_bytes(b);
    if !(1..=4096).contains(&n) {
        return Err(Error::Format(format!("bad string length {n}")));
    }
    let mut text = vec![0; n as usize];
    c.read_exact(&mut text)?;
    String::from_utf8(text).map_err(|_| Error::Format("non-UTF8 header string".into()))
}

pub fn parse_header(bytes: &[u8]) -> Result<Header> {
    let mut c = Cursor::new(bytes);
    if string(&mut c)? != "HEADER_START" {
        return Err(Error::Format("missing HEADER_START".into()));
    }
    let mut fields = BTreeMap::new();
    loop {
        let key = string(&mut c)?;
        if key == "HEADER_END" {
            break;
        }
        let value = match key.as_str() {
            "rawdatafile" | "source_name" => Value::Text(string(&mut c)?),
            "machine_id" | "telescope_id" | "data_type" | "barycentric" | "pulsarcentric"
            | "nbits" | "nsamples" | "nchans" | "nifs" | "nbeams" | "ibeam" | "signed" => {
                let mut b = [0; 4];
                c.read_exact(&mut b)?;
                Value::Int(i32::from_le_bytes(b))
            }
            "az_start" | "za_start" | "src_raj" | "src_dej" | "tstart" | "tsamp" | "fch1"
            | "foff" | "refdm" | "period" => {
                let mut b = [0; 8];
                c.read_exact(&mut b)?;
                Value::Float(f64::from_le_bytes(b))
            }
            _ => {
                return Err(Error::Format(format!(
                    "unknown keyword {key}; cannot infer its type safely"
                )));
            }
        };
        if fields.insert(key.clone(), value).is_some() {
            return Err(Error::Format(format!("duplicate {key}")));
        }
    }
    let int = |k: &str| match fields.get(k) {
        Some(Value::Int(v)) => Ok(*v),
        _ => Err(Error::Format(format!("missing integer {k}"))),
    };
    let float = |k: &str| match fields.get(k) {
        Some(Value::Float(v)) if v.is_finite() => Ok(*v),
        _ => Err(Error::Format(format!("missing/invalid float {k}"))),
    };
    let nchans = int("nchans")?;
    let nbits = int("nbits")?;
    let tsamp = float("tsamp")?;
    let fch1 = float("fch1")?;
    let foff = float("foff")?;
    if nchans <= 0 || tsamp <= 0.0 || ![1, 2, 4, 8, 16, 32].contains(&nbits) {
        return Err(Error::Format("invalid nchans/tsamp/nbits".into()));
    }
    if fields.get("nifs").is_some_and(|v| *v != Value::Int(1)) {
        return Err(Error::Format("only nifs=1 supported".into()));
    }
    if fields.get("signed").is_some_and(|v| *v != Value::Int(0)) {
        return Err(Error::Format("signed integer samples unsupported".into()));
    }
    if fields.get("data_type").is_some_and(|v| *v != Value::Int(1)) {
        return Err(Error::Format(
            "only filterbank data_type=1 supported".into(),
        ));
    }
    let h = Header {
        nchans: nchans as usize,
        nbits: nbits as u32,
        tsamp,
        fch1,
        foff,
        data_offset: c.position() as usize,
        fields,
    };
    if h.frequencies().iter().any(|f| !f.is_finite() || *f <= 0.0) {
        return Err(Error::Format("nonpositive frequency".into()));
    }
    Ok(h)
}

pub fn decode(bytes: &[u8], limit: Option<usize>) -> Result<Filterbank> {
    let header = parse_header(bytes)?;
    let payload = &bytes[header.data_offset..];
    let bits_per_frame = header
        .nchans
        .checked_mul(header.nbits as usize)
        .ok_or_else(|| Error::Format("frame size overflow".into()))?;
    let bits = payload
        .len()
        .checked_mul(8)
        .ok_or_else(|| Error::Format("payload overflow".into()))?;
    if bits % bits_per_frame != 0 {
        return Err(Error::Format(
            "truncated frame or packed padding unsupported".into(),
        ));
    }
    let total = bits / bits_per_frame;
    let nsamples = limit.unwrap_or(total).min(total);
    if nsamples < 2 {
        return Err(Error::Format("need at least two samples".into()));
    }
    let mut data = vec![0.0; nsamples * header.nchans];
    for t in 0..nsamples {
        for ch in 0..header.nchans {
            let i = t * header.nchans + ch;
            let v = match header.nbits {
                1 | 2 | 4 => {
                    // SIGPROC packs low-bit samples in least-significant-bit order.
                    let bit = i * header.nbits as usize;
                    ((payload[bit / 8] >> (bit % 8)) & ((1 << header.nbits) - 1)) as f32
                }
                8 => payload[i] as f32,
                16 => u16::from_le_bytes(payload[i * 2..i * 2 + 2].try_into().unwrap()) as f32,
                32 => f32::from_le_bytes(payload[i * 4..i * 4 + 4].try_into().unwrap()),
                _ => unreachable!(),
            };
            if !v.is_finite() {
                return Err(Error::Format(format!("nonfinite sample at {t}/{ch}")));
            }
            data[ch * nsamples + t] = v;
        }
    }
    Ok(Filterbank {
        header,
        nsamples,
        data,
    })
}

pub fn read(path: impl AsRef<Path>, limit: Option<usize>) -> Result<Filterbank> {
    let file = File::open(path)?;
    // SAFETY: read-only mapping; callers must not mutate/truncate the file while reading.
    let map = unsafe { Mmap::map(&file)? };
    decode(&map, limit)
}
pub fn read_header(path: impl AsRef<Path>) -> Result<Header> {
    let file = File::open(path)?;
    // SAFETY: read-only mapping, same immutable-input contract as read().
    let map = unsafe { Mmap::map(&file)? };
    parse_header(&map)
}
fn put_string(w: &mut impl Write, text: &str) -> Result<()> {
    w.write_all(&(text.len() as i32).to_le_bytes())?;
    w.write_all(text.as_bytes())?;
    Ok(())
}
/// Write 32-bit float filterbanks, independent of source bit depth.
pub fn write(mut w: impl Write, fb: &Filterbank) -> Result<()> {
    if fb.data.len() != fb.nsamples * fb.header.nchans {
        return Err(Error::Format("array shape mismatch".into()));
    }
    put_string(&mut w, "HEADER_START")?;
    let mut fields = fb.header.fields.clone();
    for (k, v) in [
        ("nchans", fb.header.nchans as i32),
        ("nbits", 32),
        ("nifs", 1),
        ("data_type", 1),
        ("signed", 0),
    ] {
        fields.insert(k.into(), Value::Int(v));
    }
    fields.remove("nsamples");
    for (k, v) in [
        ("tsamp", fb.header.tsamp),
        ("fch1", fb.header.fch1),
        ("foff", fb.header.foff),
    ] {
        fields.insert(k.into(), Value::Float(v));
    }
    for (k, v) in fields {
        put_string(&mut w, &k)?;
        match v {
            Value::Int(v) => w.write_all(&v.to_le_bytes())?,
            Value::Float(v) => w.write_all(&v.to_le_bytes())?,
            Value::Text(v) => put_string(&mut w, &v)?,
        }
    }
    put_string(&mut w, "HEADER_END")?;
    for t in 0..fb.nsamples {
        for c in 0..fb.header.nchans {
            w.write_all(&fb.channel(c)[t].to_le_bytes())?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn roundtrip_preserves_channel_order_and_samples() {
        let fb = Filterbank {
            header: Header {
                fields: BTreeMap::new(),
                nchans: 2,
                nbits: 32,
                tsamp: 0.001,
                fch1: 1400.0,
                foff: -1.0,
                data_offset: 0,
            },
            nsamples: 3,
            data: vec![1., 2., 3., 4., 5., 6.],
        };
        let mut bytes = Vec::new();
        write(&mut bytes, &fb).unwrap();
        let got = decode(&bytes, None).unwrap();
        assert_eq!(got.data, fb.data);
        assert_eq!(got.header.frequencies(), vec![1400., 1399.]);
        assert!(decode(&bytes[..bytes.len() - 1], None).is_err());
    }
    #[test]
    fn rejects_unknown_keyword() {
        let mut b = vec![];
        put_string(&mut b, "HEADER_START").unwrap();
        put_string(&mut b, "surprise").unwrap();
        assert!(
            parse_header(&b)
                .unwrap_err()
                .to_string()
                .contains("unknown keyword")
        );
    }
    // Literal payload bytes exercise the decoder independently of our float writer.
    fn fixture(bits: i32, payload: &[u8]) -> Vec<u8> {
        let mut bytes = vec![];
        put_string(&mut bytes, "HEADER_START").unwrap();
        for (key, value) in [
            ("nchans", 2),
            ("nbits", bits),
            ("nifs", 1),
            ("data_type", 1),
        ] {
            put_string(&mut bytes, key).unwrap();
            bytes.extend(value.to_le_bytes());
        }
        for (key, value) in [("tsamp", 0.001f64), ("fch1", 1400.), ("foff", -1.)] {
            put_string(&mut bytes, key).unwrap();
            bytes.extend(value.to_le_bytes());
        }
        put_string(&mut bytes, "HEADER_END").unwrap();
        bytes.extend(payload);
        bytes
    }
    #[test]
    fn independent_bit_depth_fixtures() {
        let cases = [
            (1, vec![0b10100110], vec![0., 1., 0., 0., 1., 0., 1., 1.]),
            (2, vec![0b11100100], vec![0., 2., 1., 3.]),
            (4, vec![0x21, 0x43], vec![1., 3., 2., 4.]),
            (8, vec![1, 2, 3, 4], vec![1., 3., 2., 4.]),
            (16, vec![1, 0, 2, 1, 3, 0, 4, 2], vec![1., 3., 258., 516.]),
            (
                32,
                vec![0, 0, 128, 63, 0, 0, 0, 64, 0, 0, 64, 64, 0, 0, 128, 64],
                vec![1., 3., 2., 4.],
            ),
        ];
        for (bits, payload, expected) in cases {
            let bytes = fixture(bits, &payload);
            assert_eq!(decode(&bytes, None).unwrap().data, expected, "nbits {bits}");
            assert_eq!(decode(&bytes, Some(2)).unwrap().nsamples, 2);
        }
        assert!(
            decode(
                &fixture(32, &[0, 0, 128, 127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                None
            )
            .is_err()
        );
        assert!(decode(&fixture(16, &[1, 0, 2]), None).is_err());
    }
}
