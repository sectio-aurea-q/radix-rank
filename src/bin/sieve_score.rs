//! sieve-score — Bewertung vieler Kandidatenkurven mit dem Radix-Kern bei wenigen Primzahlen.
//!
//! Für jede Kurve und jedes Fenster L der Liste:
//!   b_std, b_rad   Schranken mit fester Fejér-Testfunktion und mit dem Radix-Kern (unter GRH),
//!   a_std, a_rad   Ableitung der Schranke nach log N (Leitertermgewicht der jeweiligen Testfunktion),
//!   r_std, r_rad   b − a·log N: der vom Führer bereinigte Anteil (Rangsignal der Primzahlseite).
//!
//! Aufruf: sieve-score <curves.jsonl> <out.jsonl> <L1,L2,...> [K=10]
//!   curves.jsonl: Zeilen wie für radix-rank (label, logN, bad, rank, a, primes, ap); alle Kurven
//!   mit derselben Primzahlliste.
#[path = "../quad.rs"]
#[allow(dead_code)]
mod quad;
#[path = "../window.rs"]
#[allow(dead_code)]
mod window;
#[path = "../curve.rs"]
#[allow(dead_code)]
mod curve;

use curve::{prime_powers, CurveData};
use serde::Serialize;
use std::io::{BufRead, BufReader, Write};
use window::{cholesky_solve, Window};

#[derive(Serialize)]
struct Score { l: f64, b_std: f64, b_rad: f64, a_std: f64, a_rad: f64, r_std: f64, r_rad: f64 }

#[derive(Serialize)]
struct Out { label: String, #[serde(rename = "logN")] log_n: f64, rank: u32, scores: Vec<Score> }

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let ls: Vec<f64> = args[3].split(',').map(|s| s.parse().unwrap()).collect();
    let k: usize = args.get(4).map(|s| s.parse().unwrap()).unwrap_or(10);
    let rd = BufReader::new(std::fs::File::open(&args[1]).unwrap());
    let mut out = std::fs::File::create(&args[2]).unwrap();

    let mut windows: Option<(Vec<(u64, u64, u32)>, Vec<Window>)> = None;
    let mut primes_ref: Vec<u64> = Vec::new();
    let two_log_2pi = 2.0 * (2.0 * std::f64::consts::PI).ln();
    let mut n = 0usize;
    let t0 = std::time::Instant::now();
    for line in rd.lines() {
        let mut cv: CurveData = serde_json::from_str(&line.unwrap()).unwrap();
        cv.load_bin();
        if windows.is_none() {
            primes_ref = cv.primes.clone();
            let nmax = *primes_ref.last().unwrap();
            let lmax = ls.iter().cloned().fold(0.0, f64::max);
            assert!((2.0 * lmax).exp() <= nmax as f64 + 1.0, "Primzahlliste reicht nicht bis e^(2 Lmax)");
            let pp = prime_powers(&primes_ref, nmax);
            let logs: Vec<f64> = pp.iter().map(|&(m, _, _)| (m as f64).ln()).collect();
            let ws = ls.iter().map(|&l| Window::new(l, k, &logs)).collect();
            windows = Some((pp, ws));
        }
        assert_eq!(cv.primes, primes_ref, "alle Kurven brauchen dieselbe Primzahlliste");
        let (pp, ws) = windows.as_ref().unwrap();
        let w = cv.weights(pp);
        let log_n = cv.log_conductor();
        let cond = log_n - two_log_2pi;
        let kk = k * k;
        let mut scores = Vec::new();
        for win in ws {
            let mut wm: Vec<f64> = (0..kk).map(|ij| win.g0[ij] * cond + win.arch[ij]).collect();
            for (wn, gn) in w.iter().zip(&win.gn) {
                for ij in 0..kk { wm[ij] -= 2.0 * wn * gn[ij]; }
            }
            let v0 = 2.0 * win.l;
            let mut v = vec![0.0; k];
            v[0] = v0;
            let b_std = wm[0] / (v0 * v0);
            let a_std = win.g0[0] / (v0 * v0);
            let (b_rad, a_rad) = match cholesky_solve(&wm, &v, k) {
                Some(x) => {
                    let kern0: f64 = x.iter().zip(&v).map(|(a, b)| a * b).sum();
                    let c: Vec<f64> = x.iter().map(|a| a / kern0).collect();
                    let mut a = 0.0;
                    for i in 0..k { for j in 0..k { a += c[i] * win.g0[i * k + j] * c[j]; } }
                    (1.0 / kern0, a)
                }
                None => (f64::NAN, f64::NAN),
            };
            scores.push(Score { l: win.l, b_std, b_rad, a_std, a_rad,
                                r_std: b_std - a_std * log_n, r_rad: b_rad - a_rad * log_n });
        }
        writeln!(out, "{}", serde_json::to_string(&Out { label: cv.label.clone(), log_n, rank: cv.rank, scores }).unwrap()).unwrap();
        n += 1;
    }
    eprintln!("{} Kurven in {:.2?}", n, t0.elapsed());
}
