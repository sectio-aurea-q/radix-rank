//! Elliptische Kurve: Gewichte Lambda(n) c_n / sqrt(n) der logarithmischen Ableitung von L(s).
use serde::Deserialize;

#[derive(Deserialize)]
pub struct CurveData {
    pub label: String,
    #[serde(rename = "N")]
    pub n: u64,
    pub rank: u32,
    #[allow(dead_code)]
    pub a: serde_json::Value, // Weierstrass-Koeffizienten (können i64 sprengen)
    pub primes: Vec<u64>,
    pub ap: Vec<i64>,
}

/// Alle Primzahlpotenzen n = p^k <= nmax, aufsteigend: (n, p, k).
pub fn prime_powers(primes: &[u64], nmax: u64) -> Vec<(u64, u64, u32)> {
    let mut v = Vec::new();
    for &p in primes {
        let (mut q, mut k) = (p, 1u32);
        while q <= nmax {
            v.push((q, p, k));
            q = match q.checked_mul(p) { Some(x) => x, None => break };
            k += 1;
        }
    }
    v.sort();
    v
}

impl CurveData {
    /// Gewichte in der Reihenfolge von `pp` (aus prime_powers mit denselben Primzahlen).
    pub fn weights(&self, pp: &[(u64, u64, u32)]) -> Vec<f64> {
        let idx: std::collections::HashMap<u64, i64> =
            self.primes.iter().cloned().zip(self.ap.iter().cloned()).collect();
        pp.iter()
            .map(|&(n, p, k)| {
                let s1 = idx[&p] as f64 / (p as f64).sqrt();
                let ck = if self.n % p == 0 {
                    s1.powi(k as i32) // schlechte Reduktion: a_p in {0, ±1}
                } else {
                    // alpha^k + beta^k, alpha + beta = s1, alpha beta = 1
                    let (mut s_prev, mut s_cur) = (2.0, s1);
                    for _ in 1..k {
                        let s_next = s1 * s_cur - s_prev;
                        s_prev = s_cur;
                        s_cur = s_next;
                    }
                    s_cur
                };
                (p as f64).ln() * ck / (n as f64).sqrt()
            })
            .collect()
    }
}
