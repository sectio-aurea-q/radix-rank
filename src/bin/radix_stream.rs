//! radix-stream — radix-rank for very long prime lists, streamed block-wise from compact .rap files (written by `ap`).
//!
//! Nothing is held in memory but one block per thread and the a_p of the primes below 10^6 (for prime powers).
//! For a target T it finds, by bisection on the L-grid (step dL), the first window where the bound drops below T,
//! for the Radix kernel and for the fixed Fejér kernel (Bober), and also reports both bounds at the largest
//! window the data allow. Output: one JSON line; every hit carries the test-function coefficients for the Arb
//! certification (scripts/cert_rap.py).
//!
//! Usage: radix-stream <curve.json> <rap_prefix> <target> [K=10] [dL=0.01] [Lmin=5]
//!   curve.json: {"label", "logN", "bad": [primes < 2^63], "rank"}
//!   rap_prefix: e.g. data/c4239  ->  data/c4239.<lo>-<hi>.rap (+ .ok), contiguous from 0
#[path = "../quad.rs"]
#[allow(dead_code)]
mod quad;
#[path = "../window.rs"]
#[allow(dead_code)]
mod window;

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use window::{cholesky_solve, Window};

#[derive(Deserialize)]
struct Curve { label: String, #[serde(rename = "logN")] log_n: f64, bad: Vec<u64>, rank: u32 }

#[derive(Serialize, Clone)]
struct Hit { l: f64, primes: u64, bound: f64, coeffs: Vec<f64> }

#[derive(Serialize)]
struct Out {
    label: String,
    #[serde(rename = "logN")]
    log_n: f64,
    rank: u32,
    target: f64,
    pmax: u64,
    radix_hit: Option<Hit>,
    std_hit: Option<Hit>,
    radix_end: Hit,
    std_end: Hit,
    evals: Vec<(f64, u64, f64, f64)>, // (L, primes, Bober, Radix)
}

struct Block { lo: u64, hi: u64, path: String }

/// read a .rap block: (primes, a_p)
fn read_rap(path: &str) -> (Vec<u64>, Vec<i32>) {
    let v = std::fs::read(path).unwrap();
    assert_eq!(&v[0..4], b"RAP1", "{path}: not a RAP1 file");
    let rd = |o: usize| u64::from_le_bytes(v[o..o + 8].try_into().unwrap());
    let (lo, n) = (rd(4), rd(20) as usize);
    let mut ps = Vec::with_capacity(n);
    let mut p = lo;
    for i in 0..n {
        p += u16::from_le_bytes(v[28 + 2 * i..30 + 2 * i].try_into().unwrap()) as u64;
        ps.push(p);
    }
    let o = 28 + 2 * n;
    let aps = (0..n).map(|i| i32::from_le_bytes(v[o + 4 * i..o + 4 * i + 4].try_into().unwrap())).collect();
    (ps, aps)
}

fn blocks(prefix: &str) -> Vec<Block> {
    let dir = std::path::Path::new(prefix).parent().map(|d| d.to_path_buf()).unwrap_or_else(|| ".".into());
    let dir = if dir.as_os_str().is_empty() { std::path::PathBuf::from(".") } else { dir };
    let stem = std::path::Path::new(prefix).file_name().unwrap().to_string_lossy().to_string();
    let mut bs: Vec<Block> = std::fs::read_dir(&dir).unwrap().filter_map(|e| {
        let name = e.ok()?.file_name().to_string_lossy().to_string();
        let rest = name.strip_prefix(&format!("{stem}."))?.strip_suffix(".rap")?;
        let (a, b) = rest.split_once('-')?;
        if !dir.join(format!("{name}").replace(".rap", ".ok")).exists() { return None; }
        Some(Block { lo: a.parse().ok()?, hi: b.parse().ok()?, path: dir.join(&name).to_string_lossy().to_string() })
    }).collect();
    bs.sort_by_key(|b| b.lo);
    // keep the contiguous prefix starting at 0
    let mut end = 0u64;
    let mut out = Vec::new();
    for b in bs { if b.lo != end { break; } end = b.hi; out.push(b); }
    out
}

/// Weil-form pieces for window L: returns (Bober bound, Radix bound, Radix coefficients, primes used)
fn eval(l: f64, k: usize, cv: &Curve, bs: &[Block], small: &HashMap<u64, i64>, threads: usize) -> (f64, f64, Vec<f64>, u64) {
    let win = Window::new(l, k, &[]);
    let m = win.coef[0].len();
    let two_l = 2.0 * l;
    let pmax_f = two_l.exp();
    // ---- prime sum as Chebyshev moments: mom_j = sum_n w_n T_j(log n / L - 1)
    let add = |mom: &mut [f64], t: &mut [f64], w: f64, u: f64| {
        let s = u / l - 1.0;
        t[0] = 1.0;
        t[1] = s;
        for j in 2..m { t[j] = 2.0 * s * t[j - 1] - t[j - 2]; }
        for j in 0..m { mom[j] += w * t[j]; }
    };
    let use_bs: Vec<&Block> = bs.iter().filter(|b| (b.lo as f64) < pmax_f).collect();
    let next = AtomicUsize::new(0);
    let total = Mutex::new((vec![0.0f64; m], 0u64));
    std::thread::scope(|sc| {
        for _ in 0..threads {
            sc.spawn(|| {
                let mut mom = vec![0.0; m];
                let mut t = vec![0.0; m];
                let mut cnt = 0u64;
                loop {
                    let i = next.fetch_add(1, Ordering::SeqCst);
                    if i >= use_bs.len() { break; }
                    let (ps, aps) = read_rap(&use_bs[i].path);
                    for (&p, &a) in ps.iter().zip(&aps) {
                        let u = (p as f64).ln();
                        if u >= two_l { break; }
                        cnt += 1;
                        add(&mut mom, &mut t, u * a as f64 / p as f64, u); // Lambda(p) c_p / sqrt p, c_p = a_p/sqrt p
                    }
                }
                let mut g = total.lock().unwrap();
                for j in 0..m { g.0[j] += mom[j]; }
                g.1 += cnt;
            });
        }
    });
    let (mut mom, cnt) = total.into_inner().unwrap();
    // prime powers p^k, k >= 2 (p < e^L < 10^6)
    let mut t = vec![0.0; m];
    let mut ps: Vec<(&u64, &i64)> = small.iter().collect();
    ps.sort();
    for (&p, &a) in ps {
        let lp = (p as f64).ln();
        if 2.0 * lp >= two_l { break; }
        let s1 = a as f64 / (p as f64).sqrt();
        let bad = cv.bad.contains(&p);
        let (mut s_prev, mut s_cur) = (2.0, s1);
        let mut kk = 1u32;
        loop {
            // advance to k = kk + 1
            let ck = if bad { s1.powi(kk as i32 + 1) } else { let nx = s1 * s_cur - s_prev; s_prev = s_cur; s_cur = nx; nx };
            kk += 1;
            let u = kk as f64 * lp;
            if u >= two_l { break; }
            add(&mut mom, &mut t, lp * ck / (p as f64).powf(kk as f64 / 2.0), u);
        }
    }
    let kk2 = k * k;
    let cond = cv.log_n - 2.0 * (2.0 * std::f64::consts::PI).ln();
    let wm: Vec<f64> = (0..kk2).map(|ij| {
        let pp: f64 = win.coef[ij].iter().zip(&mom).map(|(a, b)| a * b).sum();
        win.g0[ij] * cond + win.arch[ij] - 2.0 * pp
    }).collect();
    let v0 = 2.0 * l;
    let mut v = vec![0.0; k];
    v[0] = v0;
    let b_std = wm[0] / (v0 * v0);
    let (b_rad, coeffs) = match cholesky_solve(&wm, &v, k) {
        Some(x) => {
            let kern0: f64 = x.iter().zip(&v).map(|(a, b)| a * b).sum();
            (1.0 / kern0, x.iter().map(|a| a / kern0).collect())
        }
        None => (f64::NAN, vec![]),
    };
    (b_std, b_rad, coeffs, cnt)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 4 {
        eprintln!("usage: radix-stream <curve.json> <rap_prefix> <target> [K=10] [dL=0.01] [Lmin=5]");
        std::process::exit(2);
    }
    let cv: Curve = serde_json::from_str(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
    let bs = blocks(&a[2]);
    assert!(!bs.is_empty(), "no contiguous .rap blocks from 0 found for {}", a[2]);
    let target: f64 = a[3].parse().unwrap();
    let k: usize = a.get(4).map(|s| s.parse().unwrap()).unwrap_or(10);
    let dl: f64 = a.get(5).map(|s| s.parse().unwrap()).unwrap_or(0.01);
    let lmin: f64 = a.get(6).map(|s| s.parse().unwrap()).unwrap_or(5.0);
    let threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    let pmax = bs.last().unwrap().hi;
    // small primes (< 10^6) with a_p, for the prime powers
    let mut small = HashMap::new();
    for b in bs.iter().filter(|b| b.lo < 1_000_000) {
        let (ps, aps) = read_rap(&b.path);
        for (&p, &x) in ps.iter().zip(&aps) { if p < 1_000_000 { small.insert(p, x as i64); } }
    }
    let gi = |l: f64| (l / dl + 1e-9).floor() as i64;
    let i_max = gi((pmax as f64).ln() / 2.0);
    let i_min = gi(lmin);
    let t0 = std::time::Instant::now();
    let cache: Arc<Mutex<HashMap<i64, (f64, f64, Vec<f64>, u64)>>> = Arc::new(Mutex::new(HashMap::new()));
    let ev = |i: i64| -> (f64, f64, Vec<f64>, u64) {
        if let Some(r) = cache.lock().unwrap().get(&i) { return r.clone(); }
        let l = (i as f64 * dl * 1e6).round() / 1e6;
        let r = eval(l, k, &cv, &bs, &small, threads);
        eprintln!("L = {l:.2}: primes {:>11}  Bober {:.6}  Radix {:.6}  ({:.1?})", r.3, r.0, r.1, t0.elapsed());
        cache.lock().unwrap().insert(i, r.clone());
        r
    };
    let l_of = |i: i64| (i as f64 * dl * 1e6).round() / 1e6;
    let end = ev(i_max);
    // first grid index in [i_min, i_max] with bound < target, by bisection (bounds decrease with L)
    let find = |which: usize| -> Option<i64> {
        let val = |r: &(f64, f64, Vec<f64>, u64)| if which == 0 { r.0 } else { r.1 };
        if !(val(&ev(i_max)) < target) { return None; }
        if val(&ev(i_min)) < target { return Some(i_min); }
        let (mut lo, mut hi) = (i_min, i_max); // val(lo) >= target, val(hi) < target
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if val(&ev(mid)) < target { hi = mid } else { lo = mid }
        }
        Some(hi)
    };
    let hit = |i: Option<i64>, which: usize| i.map(|i| {
        let r = ev(i);
        let mut e0 = vec![0.0; k];
        e0[0] = 1.0 / (2.0 * l_of(i));
        if which == 0 { Hit { l: l_of(i), primes: r.3, bound: r.0, coeffs: e0 } }
        else { Hit { l: l_of(i), primes: r.3, bound: r.1, coeffs: r.2.clone() } }
    });
    let ir = find(1);
    let is = find(0);
    let mut e0 = vec![0.0; k];
    e0[0] = 1.0 / (2.0 * l_of(i_max));
    let mut evals: Vec<(f64, u64, f64, f64)> = cache.lock().unwrap().iter().map(|(&i, r)| (l_of(i), r.3, r.0, r.1)).collect();
    evals.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let out = Out {
        label: cv.label.clone(), log_n: cv.log_n, rank: cv.rank, target, pmax,
        radix_hit: hit(ir, 1), std_hit: hit(is, 0),
        radix_end: Hit { l: l_of(i_max), primes: end.3, bound: end.1, coeffs: end.2.clone() },
        std_end: Hit { l: l_of(i_max), primes: end.3, bound: end.0, coeffs: e0 },
        evals,
    };
    println!("{}", serde_json::to_string(&out).unwrap());
}
