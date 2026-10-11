//! radix-rank — implementation of the RADIX RIEMANNIANA (S. Hackfort, 2026).
//!
//!   ℛ_{E,L}(w,z) = Φ(w)* G⁻¹ Φ(z)      reproducing kernel of the Weil form of E on the window [-L, L],
//!                                        built from Gamma factor, conductor and primes p < e^{2L} only.
//!   r_an(E) ≤ 1 / ℛ_{E,L}(0,0)          (under GRH) — the optimal windowed rank bound.
//!
//! radix-rank: Obergrenzen für den analytischen Rang elliptischer Kurven (unter GRH)
//! über den Radix-Kern (reproduzierender Kern der Weil-Form im Fenster [-L, L]).
//!
//! Aufruf: radix-rank <curves.jsonl> <out.jsonl> [Lmin Lmax dL K]
//! Gleitkomma-Phase: findet je Kurve das kleinste L, ab dem die Schranke den Rang festlegt,
//! und gibt die Koeffizienten der Testfunktion für die beweisfeste Prüfung (Arb) aus.
mod curve;
mod quad;
mod window;

use curve::{prime_powers, CurveData};
use serde::Serialize;
use std::io::{BufRead, BufReader, Write};
use window::{cholesky_solve, Window};

#[derive(Serialize, Clone)]
struct Hit {
    l: f64,
    primes: usize,
    bound: f64,
    coeffs: Vec<f64>, // Testfunktion f = sum c_k P_2k(x/L), normiert auf f^(0) = 1
}

#[derive(Serialize)]
struct Result {
    label: String,
    #[serde(rename = "N")]
    n: u64,
    rank: u32,
    std_r1: Option<Hit>,   // Standard-Testfunktion, Schranke < r+1
    radix_r1: Option<Hit>, // Radix-Kern, Schranke < r+1
    std_r2: Option<Hit>,   // mit Paritätsargument: Schranke < r+2 genügt
    radix_r2: Option<Hit>,
    curve_std: Vec<f64>,   // Schranken über das ganze L-Gitter
    curve_radix: Vec<f64>,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let inp = &args[1];
    let outp = &args[2];
    let num = |i: usize, d: f64| args.get(i).map(|s| s.parse().unwrap()).unwrap_or(d);
    let (lmin, lmax, dl) = (num(3, 0.3), num(4, 4.5), num(5, 0.05));
    let k = num(6, 10.0) as usize;

    let curves: Vec<CurveData> = BufReader::new(std::fs::File::open(inp).unwrap())
        .lines()
        .map(|l| serde_json::from_str(&l.unwrap()).unwrap())
        .collect();
    let primes = curves[0].primes.clone();
    let nmax = *primes.last().unwrap();
    assert!((2.0 * lmax).exp() <= nmax as f64 + 1.0, "Primzahlliste reicht nicht bis e^(2 Lmax)");
    let pp = prime_powers(&primes, nmax);
    let pp_logs: Vec<f64> = pp.iter().map(|&(n, _, _)| (n as f64).ln()).collect();

    let nl = ((lmax - lmin) / dl).round() as usize + 1;
    let t0 = std::time::Instant::now();
    let stream = std::env::var("RADIX_STREAM").is_ok();
    let empty: Vec<f64> = Vec::new();
    let windows: Vec<Window> = (0..nl)
        .map(|i| Window::new(((lmin + dl * i as f64) * 1e6).round() / 1e6, k, if stream { &empty } else { &pp_logs }))
        .collect();
    eprintln!("{} Fenster vorberechnet in {:.2?}", nl, t0.elapsed());
    let prime_count = |l: f64| primes.iter().take_while(|&&p| (p as f64).ln() < 2.0 * l).count();

    let mut out = std::fs::File::create(outp).unwrap();
    let t1 = std::time::Instant::now();
    for cv in &curves {
        assert_eq!(cv.primes.len(), primes.len());
        let w = cv.weights(&pp);
        let cond = (cv.n as f64).ln() - 2.0 * (2.0 * std::f64::consts::PI).ln();
        let r = cv.rank as f64;
        let (mut s1, mut x1, mut s2, mut x2) = (None, None, None, None);
        let (mut cs, mut cr) = (Vec::new(), Vec::new());
        for win in &windows {
            let kk = k * k;
            let mut wm: Vec<f64> = (0..kk).map(|ij| win.g0[ij] * cond + win.arch[ij]).collect();
            if stream {
                let pp_part = win.prime_part(&w, &pp_logs);
                for ij in 0..kk { wm[ij] -= 2.0 * pp_part[ij]; }
            } else {
                for (wn, gn) in w.iter().zip(&win.gn) {
                    for ij in 0..kk {
                        wm[ij] -= 2.0 * wn * gn[ij];
                    }
                }
            }
            let v0 = 2.0 * win.l;
            let mut v = vec![0.0; k];
            v[0] = v0;
            let b_std = wm[0] / (v0 * v0);
            let (b_rad, coeffs) = match cholesky_solve(&wm, &v, k) {
                Some(x) => {
                    let kern0: f64 = x.iter().zip(&v).map(|(a, b)| a * b).sum();
                    let c: Vec<f64> = x.iter().map(|a| a / kern0).collect(); // f^(0) = v.c = 1
                    (1.0 / kern0, c)
                }
                None => (f64::NAN, vec![]),
            };
            cs.push(b_std);
            cr.push(b_rad);
            let hit = |b: f64, c: Vec<f64>| Some(Hit { l: win.l, primes: prime_count(win.l), bound: b, coeffs: c });
            let mut e0 = vec![0.0; k];
            e0[0] = 1.0 / v0;
            if s1.is_none() && b_std < r + 1.0 { s1 = hit(b_std, e0.clone()); }
            if s2.is_none() && b_std < r + 2.0 { s2 = hit(b_std, e0.clone()); }
            if x1.is_none() && b_rad < r + 1.0 { x1 = hit(b_rad, coeffs.clone()); }
            if x2.is_none() && b_rad < r + 2.0 { x2 = hit(b_rad, coeffs.clone()); }
        }
        let res = Result { label: cv.label.clone(), n: cv.n, rank: cv.rank, std_r1: s1, radix_r1: x1,
                           std_r2: s2, radix_r2: x2, curve_std: cs, curve_radix: cr };
        writeln!(out, "{}", serde_json::to_string(&res).unwrap()).unwrap();
    }
    eprintln!("{} Kurven in {:.2?}", curves.len(), t1.elapsed());
}
