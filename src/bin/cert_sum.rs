//! cert_sum — rigorous enclosure of the prime sum of the explicit formula, streamed from .rap blocks, all cores.
//!
//!   S = sum_{p^k < e^{2L}} log p * c_{p^k} * p^{-k/2} * g(k log p),   c_{p^k} = alpha^k + beta^k (good p), (a_p/sqrt p)^k (bad p)
//!
//! g is given by its Chebyshev coefficients on [0, 2L] as f64 intervals (computed exactly from the rational
//! test function by scripts/cert_fast.py). All arithmetic is interval arithmetic in f64: every +, -, *, /, sqrt
//! is computed with IEEE round-to-nearest (exact to half an ulp) and the result is widened outward by one ulp,
//! so the true value is always contained. log p is computed from the exact reduction p = m 2^e and the atanh
//! series with a rigorous tail bound; log 2 is an interval constant checked against Arb by cert_fast.py.
//! Partial sums are kept per block (hierarchical summation keeps the widening small).
//!
//! Usage: cert_sum <job.json> <rap_prefix> [threads]     -> prints JSON {lo, hi (hex f64), primes, terms}
//!   job.json: {"cheb": [[lo,hi],...] (hex), "inv_l": [lo,hi], "two_l": [lo,hi], "ln2": [lo,hi], "bad": [...]}

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

#[derive(Clone, Copy, Debug)]
struct I { lo: f64, hi: f64 }

#[inline(always)]
fn dn(x: f64) -> f64 { x.next_down() }
#[inline(always)]
fn up(x: f64) -> f64 { x.next_up() }

impl I {
    fn pt(x: f64) -> I { I { lo: x, hi: x } }
    #[inline(always)]
    fn add(self, o: I) -> I { I { lo: dn(self.lo + o.lo), hi: up(self.hi + o.hi) } }
    #[inline(always)]
    fn sub(self, o: I) -> I { I { lo: dn(self.lo - o.hi), hi: up(self.hi - o.lo) } }
    #[inline(always)]
    fn mul(self, o: I) -> I {
        let a = self.lo * o.lo; let b = self.lo * o.hi; let c = self.hi * o.lo; let d = self.hi * o.hi;
        I { lo: dn(a.min(b).min(c).min(d)), hi: up(a.max(b).max(c).max(d)) }
    }
    #[inline(always)]
    fn scale(self, k: f64) -> I { // k exact, positive
        I { lo: dn(self.lo * k), hi: up(self.hi * k) }
    }
    fn recip_pos(self) -> I { // self > 0
        assert!(self.lo > 0.0);
        I { lo: dn(1.0 / self.hi), hi: up(1.0 / self.lo) }
    }
    fn sqrt_pos(self) -> I {
        assert!(self.lo >= 0.0);
        I { lo: dn(self.lo.sqrt()).max(0.0), hi: up(self.hi.sqrt()) }
    }
}

fn hex(s: &str) -> f64 { f64::from_bits(u64::from_str_radix(s.trim_start_matches("0x"), 16).unwrap()) }
fn pair(v: &serde_json::Value) -> I {
    let a = v.as_array().unwrap();
    I { lo: hex(a[0].as_str().unwrap()), hi: hex(a[1].as_str().unwrap()) }
}

/// rigorous enclosure of ln p for an integer 2 <= p < 2^53
fn ln_int(p: u64, ln2: I) -> I {
    let e = 63 - p.leading_zeros() as i32;              // 2^e <= p < 2^(e+1)
    let m = p as f64 / (2f64).powi(e);                  // exact, m in [1, 2)
    if m == 1.0 { return ln2.scale(e as f64); }
    // ln m = 2 atanh(z), z = (m-1)/(m+1) in (0, 1/3)
    let num = I::pt(m - 1.0);                           // exact (Sterbenz)
    let den = I { lo: dn(m + 1.0), hi: up(m + 1.0) };   // m+1 may round: enclose it
    let z = num.mul(den.recip_pos());
    let z2 = z.mul(z);
    // sum_{k=0}^{N-1} z^(2k+1)/(2k+1), tail <= z^(2N+1) / ((2N+1)(1 - z^2))
    let n = 18;
    let mut zp = z;
    let mut s = I::pt(0.0);
    for k in 0..n {
        s = s.add(zp.mul(I::pt((2 * k + 1) as f64).recip_pos()));
        zp = zp.mul(z2);
    }
    // tail bound with z <= z.hi < 1/3: zp.hi / ((2n+1) (1 - z2.hi))
    let tail = up(up(zp.hi / (2 * n + 1) as f64) / dn(1.0 - z2.hi));
    let s = s.add(I { lo: 0.0, hi: tail });
    ln2.scale(e as f64).add(s.scale(2.0))
}

/// g(u) for u in [0, 2L] from Chebyshev coefficients: s = u/L - 1, sum_j a_j T_j(s)
#[inline(always)]
fn g_eval(cheb: &[I], inv_l: I, u: I) -> I {
    let s = u.mul(inv_l).sub(I::pt(1.0));
    let s = I { lo: s.lo.max(-1.0), hi: s.hi.min(1.0) }; // u in [0,2L] is checked by the caller
    let mut t0 = I::pt(1.0);
    let mut t1 = s;
    let mut acc = cheb[0].add(cheb[1].mul(t1));
    let s2 = s.scale(2.0);
    for c in &cheb[2..] {
        let t2 = s2.mul(t1).sub(t0);
        let t2 = I { lo: t2.lo.max(-1.0), hi: t2.hi.min(1.0) }; // |T_j(s)| <= 1 on [-1,1]
        acc = acc.add(c.mul(t2));
        t0 = t1; t1 = t2;
    }
    acc
}

fn read_rap(path: &str) -> (Vec<u64>, Vec<i32>) {
    let v = std::fs::read(path).unwrap();
    assert_eq!(&v[0..4], b"RAP1");
    let rd = |o: usize| u64::from_le_bytes(v[o..o + 8].try_into().unwrap());
    let (lo, n) = (rd(4), rd(20) as usize);
    let mut ps = Vec::with_capacity(n);
    let mut p = lo;
    for i in 0..n { p += u16::from_le_bytes(v[28 + 2 * i..30 + 2 * i].try_into().unwrap()) as u64; ps.push(p); }
    let o = 28 + 2 * n;
    let aps = (0..n).map(|i| i32::from_le_bytes(v[o + 4 * i..o + 4 * i + 4].try_into().unwrap())).collect();
    (ps, aps)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() == 4 && a[1] == "lntest" {
        // cert_sum lntest <ln2_lo_hex> <ln2_hi_hex>: reads integers from stdin, prints "p lo_hex hi_hex"
        let ln2 = I { lo: hex(&a[2]), hi: hex(&a[3]) };
        let mut s = String::new();
        std::io::Read::read_to_string(&mut std::io::stdin(), &mut s).unwrap();
        for t in s.split_whitespace() {
            let p: u64 = t.parse().unwrap();
            let r = ln_int(p, ln2);
            println!("{p} {:016x} {:016x}", r.lo.to_bits(), r.hi.to_bits());
        }
        return;
    }
    let job: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
    let prefix = &a[2];
    let threads = a.get(3).map(|s| s.parse().unwrap())
        .unwrap_or_else(|| std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1));
    let cheb: Vec<I> = job["cheb"].as_array().unwrap().iter().map(pair).collect();
    let inv_l = pair(&job["inv_l"]);
    let two_l = pair(&job["two_l"]);
    let ln2 = pair(&job["ln2"]);
    let bad: Vec<u64> = job["bad"].as_array().unwrap().iter().map(|x| x.as_u64().unwrap()).collect();
    // blocks below e^{2L}
    let dir = std::path::Path::new(prefix).parent().unwrap();
    let stem = std::path::Path::new(prefix).file_name().unwrap().to_string_lossy().to_string();
    let pmax = two_l.hi.exp() * 1.000001 + 10.0;
    let mut blocks: Vec<(u64, u64, String)> = std::fs::read_dir(dir).unwrap().filter_map(|e| {
        let name = e.ok()?.file_name().to_string_lossy().to_string();
        let rest = name.strip_prefix(&format!("{stem}."))?.strip_suffix(".rap")?;
        let (x, y) = rest.split_once('-')?;
        let (lo, hi): (u64, u64) = (x.parse().ok()?, y.parse().ok()?);
        if (lo as f64) > pmax { return None; }
        if !dir.join(name.replace(".rap", ".ok")).exists() { return None; }
        Some((lo, hi, dir.join(&name).to_string_lossy().to_string()))
    }).collect();
    blocks.sort();
    let mut end = 0;
    for b in &blocks { assert_eq!(b.0, end, "gap in a_p data at {end}"); end = b.1; }
    assert!(end as f64 > two_l.hi.exp(), "a_p data end {end} below e^(2L)");
    let next = AtomicUsize::new(0);
    let parts: Mutex<Vec<(usize, I, u64, u64)>> = Mutex::new(Vec::new());
    std::thread::scope(|sc| {
        for _ in 0..threads {
            sc.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= blocks.len() { break; }
                let (ps, aps) = read_rap(&blocks[i].2);
                let mut sum = I::pt(0.0);
                let (mut np, mut nt) = (0u64, 0u64);
                for (&p, &ap) in ps.iter().zip(&aps) {
                    let lp = ln_int(p, ln2);
                    if lp.lo >= two_l.hi { break; }
                    assert!(lp.hi < two_l.lo, "prime {p} too close to the window edge");
                    np += 1;
                    let sq = I::pt(p as f64).sqrt_pos();
                    let isq = sq.recip_pos();
                    let s1 = I::pt(ap as f64).mul(isq);
                    let is_bad = bad.binary_search(&p).is_ok();
                    let (mut s_prev, mut s_cur) = (I::pt(2.0), s1);
                    let mut pk = s1; // (a_p/sqrt p)^k for bad p
                    let mut ipk = isq; // p^{-k/2}
                    let mut k = 1u32;
                    loop {
                        let ln = lp.scale(k as f64);
                        if ln.lo >= two_l.hi { break; }
                        assert!(ln.hi < two_l.lo, "prime power {p}^{k} too close to the window edge");
                        let ck = if is_bad { pk } else { s_cur };
                        sum = sum.add(lp.mul(ck).mul(ipk).mul(g_eval(&cheb, inv_l, ln)));
                        nt += 1;
                        if is_bad { pk = pk.mul(s1); } else { let nx = s1.mul(s_cur).sub(s_prev); s_prev = s_cur; s_cur = nx; }
                        ipk = ipk.mul(isq);
                        k += 1;
                    }
                }
                parts.lock().unwrap().push((i, sum, np, nt));
            });
        }
    });
    let mut parts = parts.into_inner().unwrap();
    parts.sort_by_key(|x| x.0);
    let (mut s, mut np, mut nt) = (I::pt(0.0), 0u64, 0u64);
    for (_, b, p, t) in parts { s = s.add(b); np += p; nt += t; }
    println!("{{\"lo\":\"{:016x}\",\"hi\":\"{:016x}\",\"primes\":{},\"terms\":{},\"lo_f\":{:e},\"hi_f\":{:e}}}",
             s.lo.to_bits(), s.hi.to_bits(), np, nt, s.lo, s.hi);
}
