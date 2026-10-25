//! kapsel-scan — Raum 1 der Kapsel: Kandidatensuche in einer elliptischen Faserung über Q(T).
//!
//! Familie (wie im ICARM-Leaderboard angegeben):
//!   y^2 + (L x + B) y = x (x + D)(x + E),   B = p E + q D + p q (L - p - q),
//!   L, p, q, D, E Polynome in T mit ganzen Koeffizienten.
//! Für jede Primzahl 3 <= p < B2 wird einmal a_p(t) für alle t in F_p tabelliert (Punktzählung),
//! danach kostet ein Parameter T = u/v nur noch eine Tabellenabfrage je Primzahl.
//!
//! Stufe 1: Mestre–Nagao-Summe S1 über p < B1 für alle T im Bereich, Bestenliste der Länge N.
//! Stufe 2: S1 und S2 über p < B2 für die Bestenliste, Ausgabe der besten K.
//!   S1(B) = sum_{p<B gut} (1 - (p-1)/#E(F_p)) log p,   S2(B) = (1/B) sum_{p<B gut} (-a_p log p).
//!
//! Aufruf: kapsel-scan <familie.json> <out.jsonl> <umax> <vmin> <vmax> [B1=1000] [B2=3000] [N=20000] [K=500] [threads]
//!   familie.json: {"name", "L": [Koeff. absteigend], "p": [...], "q": [...], "D": [...], "E": [...]}
//!   Parameterbereich: |u| <= umax, vmin <= v <= vmax, ggT(u, v) = 1.
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::io::Write;
use std::sync::Arc;

#[derive(Deserialize)]
struct Family { name: String, #[serde(rename = "L")] l: Vec<i128>, p: Vec<i128>, q: Vec<i128>,
                #[serde(rename = "D")] d: Vec<i128>, #[serde(rename = "E")] e: Vec<i128> }

const BAD: i16 = i16::MIN;

struct Table { p: u64, lnp: f64, ap: Vec<i16>, s1: Vec<f32> }

fn primes_upto(n: u64) -> Vec<u64> {
    let n = n as usize;
    let mut s = vec![true; n + 1];
    s[0] = false; if n >= 1 { s[1] = false; }
    let mut i = 2;
    while i * i <= n { if s[i] { let mut j = i * i; while j <= n { s[j] = false; j += i; } } i += 1; }
    (0..=n).filter(|&i| s[i]).map(|i| i as u64).collect()
}

fn ev(c: &[i128], t: i128, p: i128) -> i128 {
    let mut r = 0i128;
    for &a in c { r = (r * t + a.rem_euclid(p)) % p; }
    r
}

/// a_p(t) für alle t in F_p (p ungerade): a_p = -sum_x chi((Lx+B)^2 + 4x(x+D)(x+E)); BAD bei Delta = 0.
fn build(f: &Family, p: u64) -> Table {
    let pi = p as i128;
    let mut chi = vec![-1i8; p as usize];
    chi[0] = 0;
    for x in 1..p { chi[((x * x) % p) as usize] = 1; }
    let mut ap = vec![BAD; p as usize];
    for t in 0..pi {
        let (l, pp, q, d, e) = (ev(&f.l, t, pi), ev(&f.p, t, pi), ev(&f.q, t, pi), ev(&f.d, t, pi), ev(&f.e, t, pi));
        let b = (pp * e + q * d + pp * q % pi * ((l - pp - q).rem_euclid(pi))) % pi;
        // Weierstrass: a1=L, a2=D+E, a3=B, a4=DE, a6=0
        let (a1, a2, a3, a4) = (l, (d + e) % pi, b, d * e % pi);
        let b2 = (a1 * a1 + 4 * a2) % pi;
        let b4 = (2 * a4 + a1 * a3) % pi;
        let b6 = (a3 * a3) % pi;
        let b8 = ((a2 * a3 % pi * a3 - a1 * a3 % pi * a4 - a4 * a4) % pi + 3 * pi) % pi;
        let disc = (-(b2 * b2 % pi * b8 % pi) - 8 * (b4 * b4 % pi * b4 % pi) - 27 * b6 % pi * b6 % pi
                    + 9 * (b2 * b4 % pi * b6 % pi)).rem_euclid(pi);
        if disc == 0 { continue; }
        let mut s: i64 = 0;
        for x in 0..pi {
            let lin = (l * x + b) % pi;
            let val = (lin * lin + 4 * x * ((x + d) % pi) % pi * ((x + e) % pi)) % pi;
            s += chi[val as usize] as i64;
        }
        ap[t as usize] = (-s) as i16;
    }
    let lnp = (p as f64).ln();
    let s1 = ap.iter().map(|&a| if a == BAD { 0.0 } else {
        (1.0 - (p as f64 - 1.0) / (p as f64 + 1.0 - a as f64)) as f32 * lnp as f32 }).collect();
    Table { p, lnp, ap, s1 }
}

fn modinv(a: u64, p: u64) -> u64 {
    let (mut r, mut b, mut e) = (1u64, a % p, p - 2);
    while e > 0 { if e & 1 == 1 { r = r * b % p; } b = b * b % p; e >>= 1; }
    r
}

fn gcd(mut a: u64, mut b: u64) -> u64 { while b != 0 { let t = a % b; a = b; b = t; } a }

#[derive(Clone, Copy, PartialEq)]
struct Cand { s: f32, u: i64, v: i64 }
impl Eq for Cand {}
impl PartialOrd for Cand { fn partial_cmp(&self, o: &Self) -> Option<Ordering> { Some(self.cmp(o)) } }
// Min-Heap über s (die schwächste Kandidatin liegt oben)
impl Ord for Cand { fn cmp(&self, o: &Self) -> Ordering { o.s.partial_cmp(&self.s).unwrap_or(Ordering::Equal) } }

#[derive(Serialize)]
struct Out { family: String, u: i64, v: i64, s1_b1: f32, s1_b2: f64, s2_b2: f64, good_primes: usize }

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let fam: Family = serde_json::from_str(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
    let num = |i: usize, d: u64| a.get(i).map(|s| s.parse::<f64>().unwrap() as u64).unwrap_or(d);
    let (umax, vmin, vmax) = (num(3, 0) as i64, num(4, 1), num(5, 1));
    let (b1, b2) = (num(6, 1000), num(7, 3000));
    let (nkeep, kout) = (num(8, 20000) as usize, num(9, 500) as usize);
    let threads = num(10, std::thread::available_parallelism().map(|n| n.get() as u64).unwrap_or(4)) as usize;

    let t0 = std::time::Instant::now();
    let fam = Arc::new(fam);
    let ps: Vec<u64> = primes_upto(b2.max(b1)).into_iter().filter(|&p| p >= 3).collect();
    // Tabellen parallel bauen
    let mut tables: Vec<Option<Table>> = (0..ps.len()).map(|_| None).collect();
    std::thread::scope(|sc| {
        let chunks: Vec<Vec<(usize, u64)>> = (0..threads).map(|k| ps.iter().cloned().enumerate().skip(k).step_by(threads).collect()).collect();
        let hs: Vec<_> = chunks.into_iter().map(|ch| { let fam = fam.clone(); sc.spawn(move || ch.into_iter().map(|(i, p)| (i, build(&fam, p))).collect::<Vec<_>>()) }).collect();
        for h in hs { for (i, t) in h.join().unwrap() { tables[i] = Some(t); } }
    });
    let tables: Arc<Vec<Table>> = Arc::new(tables.into_iter().map(|t| t.unwrap()).collect());
    let n1 = tables.iter().take_while(|t| t.p < b1).count();
    eprintln!("[{}] {} Tabellen (p < {}) in {:.1?}; Stufe 1 mit {} Primzahlen, |u| <= {}, v in [{}, {}], {} Threads",
              fam.name, tables.len(), b2, t0.elapsed(), n1, umax, vmin, vmax, threads);

    // Stufe 1
    let t1 = std::time::Instant::now();
    let heaps: Vec<BinaryHeap<Cand>> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..threads).map(|k| { let tb = tables.clone(); sc.spawn(move || {
            let mut heap: BinaryHeap<Cand> = BinaryHeap::with_capacity(nkeep + 1);
            let mut res = vec![0u64; n1];
            let mut vinv = vec![0u64; n1];
            let mut v = vmin + k as u64;
            while v <= vmax {
                for i in 0..n1 {
                    let p = tb[i].p;
                    vinv[i] = if v % p == 0 { u64::MAX } else { modinv(v % p, p) };
                    res[i] = (-(umax as i128)).rem_euclid(p as i128) as u64;
                }
                for u in -umax..=umax {
                    if gcd(u.unsigned_abs(), v) == 1 {
                        let mut s = 0f32;
                        for i in 0..n1 {
                            if vinv[i] != u64::MAX { s += tb[i].s1[(res[i] * vinv[i] % tb[i].p) as usize]; }
                        }
                        if heap.len() < nkeep { heap.push(Cand { s, u, v: v as i64 }); }
                        else if s > heap.peek().unwrap().s { heap.pop(); heap.push(Cand { s, u, v: v as i64 }); }
                    }
                    for i in 0..n1 { res[i] += 1; if res[i] == tb[i].p { res[i] = 0; } }
                }
                v += threads as u64;
            }
            heap
        }) }).collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut all: Vec<Cand> = heaps.into_iter().flat_map(|h| h.into_vec()).collect();
    all.sort_by(|x, y| y.s.partial_cmp(&x.s).unwrap());
    all.truncate(nkeep);
    let ncand: u64 = (vmin..=vmax).map(|v| (-umax..=umax).filter(|u| gcd(u.unsigned_abs(), v) == 1).count() as u64).sum();
    eprintln!("Stufe 1: {} Parameter in {:.1?}", ncand, t1.elapsed());

    // Stufe 2
    let mut outs: Vec<Out> = all.iter().map(|c| {
        let (mut s1, mut s2, mut g) = (0f64, 0f64, 0usize);
        for t in tables.iter() {
            let p = t.p;
            if (c.v as u64) % p == 0 { continue; }
            let tt = ((c.u as i128).rem_euclid(p as i128) as u64) * modinv(c.v as u64 % p, p) % p;
            let x = t.ap[tt as usize];
            if x == BAD { continue; }
            s1 += t.s1[tt as usize] as f64; s2 += -(x as f64) * t.lnp; g += 1;
        }
        Out { family: fam.name.clone(), u: c.u, v: c.v, s1_b1: c.s, s1_b2: s1, s2_b2: s2 / b2 as f64, good_primes: g }
    }).collect();
    outs.sort_by(|x, y| y.s1_b2.partial_cmp(&x.s1_b2).unwrap());
    outs.truncate(kout);
    let mut f = std::fs::File::create(&a[2]).unwrap();
    for o in &outs { writeln!(f, "{}", serde_json::to_string(o).unwrap()).unwrap(); }
    eprintln!("Fertig in {:.1?}: beste {} nach S1(p < {}) in {}", t0.elapsed(), outs.len(), b2, a[2]);
}
