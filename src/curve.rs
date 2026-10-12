//! Elliptische Kurve: Gewichte Lambda(n) c_n / sqrt(n) der logarithmischen Ableitung von L(s).
use serde::Deserialize;

#[derive(Deserialize)]
pub struct CurveData {
    pub label: String,
    #[serde(rename = "N", default)]
    pub n: u64,
    /// log N for conductors beyond u64 (high-rank curves); overrides N when present.
    #[serde(rename = "logN", default)]
    pub log_n: Option<f64>,
    /// primes of bad reduction (needed when N does not fit into u64).
    #[serde(default)]
    pub bad: Option<Vec<u64>>,
    pub rank: u32,
    #[allow(dead_code)]
    pub a: serde_json::Value, // Weierstrass-Koeffizienten (können i64 sprengen)
    #[serde(default)]
    pub primes: Vec<u64>,
    #[serde(default)]
    pub ap: Vec<i64>,
    /// optional raw little-endian files (u64 primes, i64 a_p) for very long prime lists
    #[serde(default)]
    pub bin_p: Option<String>,
    #[serde(default)]
    pub bin_ap: Option<String>,
    /// optional cap on the primes read from the binary files
    #[serde(default)]
    pub pmax: Option<u64>,
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
    pub fn load_bin(&mut self) {
        if let (Some(fp), Some(fa)) = (&self.bin_p, &self.bin_ap) {
            let rp = std::fs::read(fp).unwrap();
            let ra = std::fs::read(fa).unwrap();
            let n = (rp.len() / 8).min(ra.len() / 8);
            let cap = self.pmax.unwrap_or(u64::MAX);
            self.primes.clear(); self.ap.clear();
            for i in 0..n {
                let p = u64::from_le_bytes(rp[8 * i..8 * i + 8].try_into().unwrap());
                if p > cap { break; }
                self.primes.push(p);
                self.ap.push(i64::from_le_bytes(ra[8 * i..8 * i + 8].try_into().unwrap()));
            }
        }
    }

    pub fn log_conductor(&self) -> f64 {
        self.log_n.unwrap_or_else(|| (self.n as f64).ln())
    }

    /// Gewichte in der Reihenfolge von `pp` (aus prime_powers mit denselben Primzahlen).
    pub fn weights(&self, pp: &[(u64, u64, u32)]) -> Vec<f64> {
        // primes are sorted: binary search instead of a hash map (memory-lean for 10^7+ primes)
        let ap_of = |p: u64| -> i64 { self.ap[self.primes.binary_search(&p).expect("prime missing")] };
        pp.iter()
            .map(|&(n, p, k)| {
                let s1 = ap_of(p) as f64 / (p as f64).sqrt();
                let is_bad = match &self.bad { Some(b) => b.contains(&p), None => self.n % p == 0 };
                let ck = if is_bad {
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
