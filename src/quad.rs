//! Quadratur und Polynome: Gauss-Legendre, Legendre-Polynome, Chebyshev-Interpolation.

/// Gauss-Legendre-Knoten und -Gewichte auf [-1, 1] (Newton-Iteration).
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    for i in 0..(n + 1) / 2 {
        let mut z = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 0.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, z);
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * z * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            dp = n as f64 * (z * p1 - p0) / (z * z - 1.0);
            let dz = p1 / dp;
            z -= dz;
            if dz.abs() < 1e-16 {
                break;
            }
        }
        x[i] = -z;
        x[n - 1 - i] = z;
        w[i] = 2.0 / ((1.0 - z * z) * dp * dp);
        w[n - 1 - i] = w[i];
    }
    (x, w)
}

/// Gerade Legendre-Polynome P_0, P_2, ..., P_{2K-2} an t (0 außerhalb [-1,1]).
pub fn legendre_even(t: f64, k_count: usize, out: &mut [f64]) {
    if t.abs() > 1.0 {
        out.iter_mut().for_each(|v| *v = 0.0);
        return;
    }
    let (mut p0, mut p1) = (1.0, t);
    out[0] = 1.0;
    let maxdeg = 2 * k_count - 2;
    for n in 2..=maxdeg {
        let p2 = ((2 * n - 1) as f64 * t * p1 - (n - 1) as f64 * p0) / n as f64;
        p0 = p1;
        p1 = p2;
        if n % 2 == 0 {
            out[n / 2] = p2;
        }
    }
}

/// Chebyshev-Punkte erster Art auf [-1, 1].
pub fn cheb_nodes(m: usize) -> Vec<f64> {
    (0..m)
        .map(|j| (std::f64::consts::PI * (j as f64 + 0.5) / m as f64).cos())
        .collect()
}

/// Chebyshev-Koeffizienten aus Werten an cheb_nodes(m) (exakt für Polynome vom Grad < m).
pub fn cheb_coeffs(vals: &[f64]) -> Vec<f64> {
    let m = vals.len();
    (0..m)
        .map(|k| {
            let s: f64 = (0..m)
                .map(|j| vals[j] * (std::f64::consts::PI * k as f64 * (j as f64 + 0.5) / m as f64).cos())
                .sum();
            if k == 0 { s / m as f64 } else { 2.0 * s / m as f64 }
        })
        .collect()
}

/// Clenshaw-Auswertung einer Chebyshev-Reihe bei x in [-1, 1].
pub fn cheb_eval(c: &[f64], x: f64) -> f64 {
    let (mut b1, mut b2) = (0.0, 0.0);
    for &ck in c.iter().skip(1).rev() {
        let b0 = 2.0 * x * b1 - b2 + ck;
        b2 = b1;
        b1 = b0;
    }
    x * b1 - b2 + c[0]
}
