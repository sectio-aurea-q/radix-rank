//! Fenster [-L, L]: kurvenunabhängige Bausteine der Weil-Form in der Basis
//! b_k(x) = P_{2k}(x/L), k = 0..K-1.
//!
//! W_E(f) = g(0)(log N - 2 log 2pi) + Arch(g) - 2 sum_n Lambda(n) c_n n^{-1/2} g(log n),
//! g = f * f,  Arch(g) = 2[-gamma g(0) + int_0^{2L} (g(0)-g(u)) / (e^u - 1) du - g(0) log(1-e^{-2L})].
use crate::quad::*;

pub const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;

pub struct Window {
    pub l: f64,
    pub k: usize,
    pub g0: Vec<f64>,        // K*K, g_ij(0)
    pub arch: Vec<f64>,      // K*K, archimedischer Term
    pub gn: Vec<Vec<f64>>,   // pro Primzahlpotenz n (log n < 2L): g_ij(log n), K*K
    pub coef: Vec<Vec<f64>>, // Chebyshev-Koeffizienten von g_ij auf [0, 2L]
}

/// g_ij(u) = int b_i(x) b_j(x-u) dx, für 0 <= u <= 2L.
fn corr_at(l: f64, k: usize, u: f64, gx: &[f64], gw: &[f64]) -> Vec<f64> {
    let (lo, hi) = (u - l, l);
    let mut g = vec![0.0; k * k];
    if hi <= lo {
        return g;
    }
    let (mid, rad) = (0.5 * (hi + lo), 0.5 * (hi - lo));
    let mut bi = vec![0.0; k];
    let mut bj = vec![0.0; k];
    for (x0, w0) in gx.iter().zip(gw) {
        let x = mid + rad * x0;
        let w = rad * w0;
        legendre_even(x / l, k, &mut bi);
        legendre_even((x - u) / l, k, &mut bj);
        for i in 0..k {
            for j in 0..k {
                g[i * k + j] += w * bi[i] * bj[j];
            }
        }
    }
    g
}

impl Window {
    /// `pp_logs`: log n aller Primzahlpotenzen (aufsteigend), kurvenunabhängig.
    pub fn new(l: f64, k: usize, pp_logs: &[f64]) -> Self {
        let (gx, gw) = gauss_legendre(4 * k + 20);
        // g_ij(u) ist auf [0, 2L] ein Polynom vom Grad <= 4K-1 -> exakte Chebyshev-Darstellung
        let m = 4 * k + 8;
        let nodes = cheb_nodes(m);
        let samples: Vec<Vec<f64>> = nodes.iter().map(|s| corr_at(l, k, l * (s + 1.0), &gx, &gw)).collect();
        let coef: Vec<Vec<f64>> = (0..k * k)
            .map(|ij| cheb_coeffs(&samples.iter().map(|v| v[ij]).collect::<Vec<_>>()))
            .collect();
        let g_at = |u: f64| -> Vec<f64> { g_eval(&coef, l, u) };
        let g0 = corr_at(l, k, 0.0, &gx, &gw);
        // archimedischer Term: Gauss-Panels auf [0, 2L]
        let (ax, aw) = gauss_legendre(40);
        let panels = 60;
        let mut integral = vec![0.0; k * k];
        for p in 0..panels {
            let a = 2.0 * l * p as f64 / panels as f64;
            let b = 2.0 * l * (p + 1) as f64 / panels as f64;
            for (x0, w0) in ax.iter().zip(&aw) {
                let u = 0.5 * (b + a) + 0.5 * (b - a) * x0;
                let w = 0.5 * (b - a) * w0 / u.exp_m1();
                let gu = g_at(u);
                for ij in 0..k * k {
                    integral[ij] += w * (g0[ij] - gu[ij]);
                }
            }
        }
        let tail = (-(-2.0 * l).exp_m1()).ln(); // log(1 - e^{-2L})
        let arch: Vec<f64> = (0..k * k)
            .map(|ij| 2.0 * (-EULER_GAMMA * g0[ij] + integral[ij] - g0[ij] * tail))
            .collect();
        let gn: Vec<Vec<f64>> = pp_logs.iter().take_while(|&&t| t < 2.0 * l).map(|&t| g_at(t)).collect();
        Window { l, k, g0, arch, gn, coef }
    }

    /// Streaming: Primzahlsumme direkt auswerten (für große Primzahllisten ohne Vorberechnung).
    pub fn prime_part(&self, weights: &[f64], logs: &[f64]) -> Vec<f64> {
        let kk = self.k * self.k;
        let mut acc = vec![0.0; kk];
        for (w, &t) in weights.iter().zip(logs) {
            if t >= 2.0 * self.l { break; }
            let g = g_eval(&self.coef, self.l, t);
            for ij in 0..kk { acc[ij] += w * g[ij]; }
        }
        acc
    }
}

fn g_eval(coef: &[Vec<f64>], l: f64, u: f64) -> Vec<f64> {
    let s = u / l - 1.0;
    coef.iter().map(|c| cheb_eval(c, s)).collect()
}

/// Cholesky-Lösung W x = v; None, wenn W nicht positiv definit ist.
pub fn cholesky_solve(w: &[f64], v: &[f64], k: usize) -> Option<Vec<f64>> {
    let mut lmat = vec![0.0; k * k];
    for i in 0..k {
        for j in 0..=i {
            let mut s = w[i * k + j];
            for m in 0..j {
                s -= lmat[i * k + m] * lmat[j * k + m];
            }
            if i == j {
                if s <= 0.0 {
                    return None;
                }
                lmat[i * k + i] = s.sqrt();
            } else {
                lmat[i * k + j] = s / lmat[j * k + j];
            }
        }
    }
    let mut y = vec![0.0; k];
    for i in 0..k {
        let mut s = v[i];
        for m in 0..i {
            s -= lmat[i * k + m] * y[m];
        }
        y[i] = s / lmat[i * k + i];
    }
    let mut x = vec![0.0; k];
    for i in (0..k).rev() {
        let mut s = y[i];
        for m in i + 1..k {
            s -= lmat[m * k + i] * x[m];
        }
        x[i] = s / lmat[i * k + i];
    }
    Some(x)
}
