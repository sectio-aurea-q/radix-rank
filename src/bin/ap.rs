//! ap — fast a_p for one elliptic curve over all primes in a range (radix-rank companion).
//!
//! Method (Shanks–Mestre): for a good prime p >= P_SMALL, pick x0, d = f(x0) != 0; the point (d x0, d^2) lies on the
//! twist E_d : y^2 = x^3 + A d^2 x + B d^3, and #E_d(F_p) = p + 1 - chi(d) a_p. Baby-step/giant-step over the WHOLE
//! Hasse interval |k| <= 2 sqrt(p) collects every k with (p+1-k) P = O; candidates are intersected over points until
//! exactly one a_p remains. Nothing is guessed: if the candidate set never becomes a singleton the program aborts.
//! Primes below P_SMALL and bad primes come from a table (PARI) given on the command line.
//!
//! Arithmetic: Montgomery form mod p (p < 2^62), Jacobian coordinates, mixed additions with a fixed affine point,
//! batched normalisation (one inversion per batch). Threads: blocks of the range are distributed over all cores.
//!
//! Usage: ap <c4> <c6> <small_ap.txt> <lo> <hi> <out_prefix> [threads]
//!   small_ap.txt: lines "p a_p" for all p < P_SMALL and every bad prime.
//!   Output (block-wise, resumable): <out_prefix>.<lo>-<hi>.rap (compact: prime gaps u16 + a_p i32, see write_rap)
//!   plus a .ok marker; AP_LEGACY=1 also writes .p.u64/.ap.i64. Block size AP_BLOCK (default 1e8).
//!   Two speed-ups over plain BSGS: a_p mod 2 from the 2-torsion test (halves the search interval) and affine
//!   arithmetic with batched inversions (baby steps in doubling rounds, giant steps in parallel chains).

use std::collections::HashMap;
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

const P_SMALL: u64 = 1000;

// ---------------------------------------------------------------- Montgomery arithmetic
#[derive(Clone, Copy)]
struct Mont {
    p: u64,
    pinv: u64, // -p^{-1} mod 2^64
    r2: u64,   // 2^128 mod p
    one: u64,
}

impl Mont {
    fn new(p: u64) -> Self {
        let mut inv: u64 = 1;
        for _ in 0..6 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
        }
        let r = ((1u128 << 64) % p as u128) as u64;
        let r2 = ((r as u128 * r as u128) % p as u128) as u64;
        Mont { p, pinv: inv.wrapping_neg(), r2, one: r }
    }
    /// Montgomery reduction, branch-free (p < 2^62, so t + m p < 2^127 never overflows).
    #[inline(always)]
    fn redc(&self, t: u128) -> u64 {
        let m = (t as u64).wrapping_mul(self.pinv);
        let r = ((t + m as u128 * self.p as u128) >> 64) as u64; // r < 2p
        Self::csub(r, self.p)
    }
    /// x < 2p -> x mod p without branches
    #[inline(always)]
    fn csub(x: u64, p: u64) -> u64 {
        let (d, borrow) = x.overflowing_sub(p);
        d.wrapping_add(p & 0u64.wrapping_sub(borrow as u64))
    }
    #[inline(always)]
    fn mul(&self, a: u64, b: u64) -> u64 { self.redc(a as u128 * b as u128) }
    #[inline(always)]
    fn add(&self, a: u64, b: u64) -> u64 { Self::csub(a + b, self.p) }
    #[inline(always)]
    fn sub(&self, a: u64, b: u64) -> u64 {
        let (d, borrow) = a.overflowing_sub(b);
        d.wrapping_add(self.p & 0u64.wrapping_sub(borrow as u64))
    }
    #[inline(always)]
    fn neg(&self, a: u64) -> u64 { if a == 0 { 0 } else { self.p - a } }
    fn to(&self, a: u64) -> u64 { self.mul(a % self.p, self.r2) }
    fn from(&self, a: u64) -> u64 { self.redc(a as u128) }
    /// b^e, branch-free in the bits of e (they change with every prime and would be mispredicted)
    fn pow(&self, mut b: u64, mut e: u64) -> u64 {
        let mut r = self.one;
        while e > 0 {
            let mask = 0u64.wrapping_sub(e & 1);
            r = self.mul(r, (b & mask) | (self.one & !mask));
            b = self.mul(b, b);
            e >>= 1;
        }
        r
    }
    fn inv(&self, a: u64) -> u64 { self.pow(a, self.p - 2) }
}

/// Jacobi symbol (a / n), n odd.
fn jacobi(mut a: u64, mut n: u64) -> i32 {
    a %= n;
    let mut t = 1;
    while a != 0 {
        let z = a.trailing_zeros();
        a >>= z;
        if z & 1 == 1 && (n % 8 == 3 || n % 8 == 5) { t = -t; }
        if a % 4 == 3 && n % 4 == 3 { t = -t; }
        std::mem::swap(&mut a, &mut n);
        a %= n;
    }
    if n == 1 { t } else { 0 }
}

/// Signed big integer as base-2^32 limbs (most significant first), for fast reduction mod p.
#[derive(Clone)]
struct Big { neg: bool, limbs: Vec<u32> }

impl Big {
    fn parse(s: &str) -> Big {
        let (neg, digits) = match s.strip_prefix('-') { Some(d) => (true, d), None => (false, s) };
        let mut limbs: Vec<u32> = vec![0]; // little-endian during construction
        for &ch in digits.as_bytes() {
            let mut carry = (ch - b'0') as u64;
            for l in limbs.iter_mut() { let v = *l as u64 * 10 + carry; *l = v as u32; carry = v >> 32; }
            if carry > 0 { limbs.push(carry as u32); }
        }
        limbs.reverse();
        Big { neg, limbs }
    }
    #[inline]
    fn modp(&self, p: u64) -> u64 {
        let mut r: u64 = 0;
        for &l in &self.limbs { r = (((r as u128) << 32 | l as u128) % p as u128) as u64; }
        if self.neg && r != 0 { p - r } else { r }
    }
}

// ---------------------------------------------------------------- curve arithmetic (Jacobian, y^2 = x^3 + a x + b)
#[derive(Clone, Copy)]
struct Jac { x: u64, y: u64, z: u64 } // z == 0 <=> point at infinity

struct Curve<'a> { m: &'a Mont, a: u64 }

impl<'a> Curve<'a> {
    fn dbl(&self, q: Jac) -> Jac {
        let m = self.m;
        if q.z == 0 || q.y == 0 { return Jac { x: m.one, y: m.one, z: 0 }; }
        let xx = m.mul(q.x, q.x);
        let yy = m.mul(q.y, q.y);
        let yyyy = m.mul(yy, yy);
        let zz = m.mul(q.z, q.z);
        let xyy = m.mul(q.x, yy);
        let s = m.add(m.add(xyy, xyy), m.add(xyy, xyy)); // 4 X YY
        let zzzz = m.mul(zz, zz);
        let mm = m.add(m.add(m.add(xx, xx), xx), m.mul(self.a, zzzz));
        let x3 = m.sub(m.mul(mm, mm), m.add(s, s));
        let y8 = { let t = m.add(yyyy, yyyy); let t = m.add(t, t); m.add(t, t) };
        let y3 = m.sub(m.mul(mm, m.sub(s, x3)), y8);
        let yz = m.mul(q.y, q.z);
        Jac { x: x3, y: y3, z: m.add(yz, yz) }
    }
    /// q + (x2, y2) with (x2, y2) affine.
    fn madd(&self, q: Jac, x2: u64, y2: u64) -> Jac {
        let m = self.m;
        if q.z == 0 { return Jac { x: x2, y: y2, z: m.one }; }
        let z1z1 = m.mul(q.z, q.z);
        let u2 = m.mul(x2, z1z1);
        let s2 = m.mul(y2, m.mul(q.z, z1z1));
        let h = m.sub(u2, q.x);
        let r = m.sub(s2, q.y);
        if h == 0 {
            if r == 0 { return self.dbl(q); }
            return Jac { x: m.one, y: m.one, z: 0 };
        }
        let hh = m.mul(h, h);
        let hhh = m.mul(h, hh);
        let v = m.mul(q.x, hh);
        let x3 = m.sub(m.sub(m.mul(r, r), hhh), m.add(v, v));
        let y3 = m.sub(m.mul(r, m.sub(v, x3)), m.mul(q.y, hhh));
        let z3 = m.mul(q.z, h);
        Jac { x: x3, y: y3, z: z3 }
    }
    fn scalar(&self, k: u64, x: u64, y: u64) -> Jac {
        let mut r = Jac { x: self.m.one, y: self.m.one, z: 0 };
        for bit in (0..64 - k.leading_zeros()).rev() {
            r = self.dbl(r);
            if (k >> bit) & 1 == 1 { r = self.madd(r, x, y); }
        }
        r
    }
    /// affine coordinates of finite points (batched inversion); infinity -> None
    fn normalize(&self, pts: &[Jac], out: &mut Vec<Option<(u64, u64)>>, scratch: &mut Vec<u64>) {
        let m = self.m;
        out.clear();
        scratch.clear();
        let mut acc = m.one;
        for q in pts {
            scratch.push(acc);
            if q.z != 0 { acc = m.mul(acc, q.z); }
        }
        let mut inv = m.inv(acc);
        out.resize(pts.len(), None);
        for i in (0..pts.len()).rev() {
            let q = pts[i];
            if q.z == 0 { continue; }
            let zi = m.mul(inv, scratch[i]);
            inv = m.mul(inv, q.z);
            let zi2 = m.mul(zi, zi);
            out[i] = Some((m.mul(q.x, zi2), m.mul(q.y, m.mul(zi2, zi))));
        }
    }
}

fn isqrt(n: u64) -> u64 {
    let mut r = (n as f64).sqrt() as u64;
    while r * r > n { r -= 1; }
    while (r + 1) * (r + 1) <= n { r += 1; }
    r
}

/// Scratch buffers reused across primes (one set per thread).
#[derive(Default)]
struct Buf { jac: Vec<Jac>, aff: Vec<Option<(u64, u64)>>, scr: Vec<u64>, keys: Vec<u64>, vals: Vec<u32>, ys: Vec<u64>,
             bx: Vec<u64>, by: Vec<u64>, sx: Vec<u64>, sy: Vec<u64>, dx: Vec<u64>, iv: Vec<u64>, pre: Vec<u64>,
             n_fast: usize, n_slow: usize }

const EMPTY: u64 = u64::MAX;
impl Buf {
    fn table_reset(&mut self, n: usize) {
        let size = (4 * n + 4).next_power_of_two();
        self.keys.clear(); self.keys.resize(size, EMPTY);
        self.vals.clear(); self.vals.resize(size, 0);
        self.ys.clear(); self.ys.resize(size, 0);
    }
    #[inline(always)]
    fn slot(&self, x: u64) -> usize { (x.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 32) as usize & (self.keys.len() - 1) }
    /// insert; returns false if x was already present
    fn insert(&mut self, x: u64, y: u64, v: u32) -> bool {
        let mut i = self.slot(x);
        loop {
            if self.keys[i] == EMPTY { self.keys[i] = x; self.vals[i] = v; self.ys[i] = y; return true; }
            if self.keys[i] == x { return false; }
            i = (i + 1) & (self.keys.len() - 1);
        }
    }
    #[inline(always)]
    fn get(&self, x: u64) -> Option<(u32, u64)> {
        let mut i = self.slot(x);
        loop {
            let k = self.keys[i];
            if k == EMPTY { return None; }
            if k == x { return Some((self.vals[i], self.ys[i])); }
            i = (i + 1) & (self.keys.len() - 1);
        }
    }
}

/// Does x^3 + a x + b have a root in F_p?  (gcd with x^p - x; a, b in Montgomery form)
fn cubic_has_root(m: &Mont, a: u64, b: u64) -> bool {
    // r = x^p mod f, f = x^3 + a x + b, as r0 + r1 x + r2 x^2
    let (na, nb) = (m.neg(a), m.neg(b));
    let (mut r0, mut r1, mut r2) = (m.one, 0u64, 0u64);
    let p = m.p;
    for bit in (0..64 - p.leading_zeros()).rev() {
        // square: d0..d4, then x^4 = -a x^2 - b x, x^3 = -a x - b
        let d0 = m.mul(r0, r0);
        let t01 = m.mul(r0, r1); let d1 = m.add(t01, t01);
        let t02 = m.mul(r0, r2); let d2 = m.add(m.mul(r1, r1), m.add(t02, t02));
        let t12 = m.mul(r1, r2); let d3 = m.add(t12, t12);
        let d4 = m.mul(r2, r2);
        r2 = m.add(d2, m.mul(na, d4));
        r1 = m.add(m.add(d1, m.mul(nb, d4)), m.mul(na, d3));
        r0 = m.add(d0, m.mul(nb, d3));
        if (p >> bit) & 1 == 1 {
            // times x: r0 x + r1 x^2 + r2 x^3 = -b r2 + (r0 - a r2) x + r1 x^2
            let (n0, n1, n2) = (m.mul(nb, r2), m.add(r0, m.mul(na, r2)), r1);
            r0 = n0; r1 = n1; r2 = n2;
        }
    }
    let (g0, g1, g2) = (r0, m.sub(r1, m.one), r2); // g = x^p - x mod f
    let f_at = |x: u64| m.add(m.add(m.mul(m.mul(x, x), x), m.mul(a, x)), b);
    if g2 != 0 {
        let iv = m.inv(g2);
        let (g1, g0) = (m.mul(g1, iv), m.mul(g0, iv)); // monic g = x^2 + g1 x + g0
        // f mod g = (g1^2 - g0 + a) x + (g1 g0 + b)
        let u1 = m.add(m.sub(m.mul(g1, g1), g0), a);
        let u0 = m.add(m.mul(g1, g0), b);
        if u1 == 0 { return u0 == 0; }
        let x = m.neg(m.mul(u0, m.inv(u1)));
        return m.add(m.add(m.mul(x, x), m.mul(g1, x)), g0) == 0;
    }
    if g1 != 0 { return f_at(m.neg(m.mul(g0, m.inv(g1)))) == 0; }
    g0 == 0
}

/// Batched inversion of nonzero values (Montgomery's trick): out[i] = vals[i]^{-1}.
#[inline]
fn batch_inv(m: &Mont, vals: &[u64], out: &mut Vec<u64>, pre: &mut Vec<u64>) {
    pre.clear();
    let mut acc = m.one;
    for &v in vals { pre.push(acc); acc = m.mul(acc, v); }
    let mut inv = m.inv(acc);
    out.clear();
    out.resize(vals.len(), 0);
    for i in (0..vals.len()).rev() {
        out[i] = m.mul(inv, pre[i]);
        inv = m.mul(inv, vals[i]);
    }
}

const CHAINS: usize = 32;
static FAST: AtomicUsize = AtomicUsize::new(0);
static SLOW: AtomicUsize = AtomicUsize::new(0);


/// Fast path of `kset`: affine arithmetic with batched inversions (baby steps in doubling rounds, giant steps in
/// CHAINS parallel chains). Returns None on ANY special case (equal x in an addition, point of small order,
/// point at infinity on the way) — the caller then uses the robust Jacobian `kset`. Same output otherwise.
fn kset_fast(c: &Curve, p: u64, b2: i64, px: u64, py: u64, e: i64, st: i64, buf: &mut Buf) -> Option<Vec<i64>> {
    // searches k = e + st * t in [-b2, b2] (st = 1, or st = 2 when the parity e of k is known) with (p+1-k) P = O,
    // as Q = t R with R = st P, Q = (p+1-e) P
    let m = c.m;
    if py == 0 { return None; }
    let dbl_aff = |x: u64, y: u64| -> Option<(u64, u64)> {
        if y == 0 { return None; }
        let num = m.add(m.add(m.add(m.mul(x, x), m.mul(x, x)), m.mul(x, x)), c.a);
        let lam = m.mul(num, m.inv(m.add(y, y)));
        let x2 = m.sub(m.sub(m.mul(lam, lam), x), x);
        Some((x2, m.sub(m.mul(lam, m.sub(x, x2)), y)))
    };
    let (rx, ry) = if st == 1 { (px, py) } else { dbl_aff(px, py)? };
    let t_lo = (-b2 - e).div_euclid(st) + if (-b2 - e).rem_euclid(st) != 0 { 1 } else { 0 };
    let t_hi = (b2 - e).div_euclid(st);
    let w = (t_hi - t_lo + 1).max(1) as u64;
    let mm = ((w as f64 / 2.0).sqrt().ceil() as usize).max(2);
    let h = 2 * mm as u64 + 1;
    // ---- baby steps B_i = (i+1) R, i = 0..mm  (one extra: B_mm = (mm+1) R, used for H)
    let nb = mm + 1;
    buf.bx.clear(); buf.by.clear();
    buf.bx.resize(nb, 0); buf.by.resize(nb, 0);
    buf.bx[0] = rx; buf.by[0] = ry;
    let (x2, y2) = dbl_aff(rx, ry)?;
    buf.bx[1] = x2;
    buf.by[1] = y2;
    let mut t = 2usize; // B_0..B_{t-1} known
    while t < nb {
        let (tx, ty) = (buf.bx[t - 1], buf.by[t - 1]); // tP; B_i + tP = B_{t+i} for i = 0..t-2 (i = t-1 would double)
        let n = (t - 1).min(nb - t);
        buf.dx.clear();
        for i in 0..n {
            let d = m.sub(tx, buf.bx[i]);
            if d == 0 { return None; }
            buf.dx.push(d);
        }
        batch_inv(m, &buf.dx, &mut buf.iv, &mut buf.pre);
        for i in 0..n {
            let l = m.mul(m.sub(ty, buf.by[i]), buf.iv[i]);
            let x3 = m.sub(m.sub(m.mul(l, l), buf.bx[i]), tx);
            buf.bx[t + i] = x3;
            buf.by[t + i] = m.sub(m.mul(l, m.sub(buf.bx[i], x3)), buf.by[i]);
        }
        t += n;
    }
    buf.table_reset(mm);
    for i in 0..mm {
        if !buf.insert(buf.bx[i], buf.by[i], i as u32 + 1) { return None; } // order <= 2 mm
    }
    // ---- H = h P = mm P + (mm+1) P
    let (ax, ay, bxm, bym) = (buf.bx[mm - 1], buf.by[mm - 1], buf.bx[mm], buf.by[mm]);
    let d = m.sub(bxm, ax);
    if d == 0 { return None; }
    let l = m.mul(m.sub(bym, ay), m.inv(d));
    let hx = m.sub(m.sub(m.mul(l, l), ax), bxm);
    let hy = m.sub(m.mul(l, m.sub(ax, hx)), ay);
    let nhy = m.neg(hy);
    // ---- giant steps G_i = Q - (c0 + i h) R = (p+1-e-st(c0+i h)) P, split into CHAINS chains of length ln
    let c0 = t_lo + mm as i64;
    let steps = ((w as u64) / h + 2) as usize;
    let nch = CHAINS.min(steps);
    let ln = (steps + nch - 1) / nch;
    // chain starts S_r = G_{r ln} = G_0 - r (ln H)
    let g0 = c.scalar((p as i64 + 1 - e - st * c0) as u64, px, py);
    let lh = c.scalar(ln as u64, hx, hy);
    buf.jac.clear();
    buf.jac.push(lh);
    c.normalize(&buf.jac, &mut buf.aff, &mut buf.scr);
    let (lhx, lhy) = buf.aff[0]?;
    let nlhy = m.neg(lhy);
    buf.jac.clear();
    let mut g = g0;
    for _ in 0..nch { buf.jac.push(g); g = c.madd(g, lhx, nlhy); }
    c.normalize(&buf.jac, &mut buf.aff, &mut buf.scr);
    buf.sx.clear(); buf.sy.clear();
    for r in 0..nch { let (x, y) = buf.aff[r]?; buf.sx.push(x); buf.sy.push(y); }
    let mut out = Vec::new();
    for s in 0..ln {
        for r in 0..nch {
            let i = r * ln + s;
            if i >= steps { continue; }
            if let Some((j, yj)) = buf.get(buf.sx[r]) {
                let base = c0 + i as i64 * h as i64;
                let t = if buf.sy[r] == yj { base + j as i64 } else { base - j as i64 };
                out.push(e + st * t);
            }
        }
        if s + 1 == ln { break; }
        // S_r <- S_r - H for all chains (batched)
        buf.dx.clear();
        for r in 0..nch {
            let d = m.sub(hx, buf.sx[r]);
            if d == 0 { return None; }
            buf.dx.push(d);
        }
        batch_inv(m, &buf.dx, &mut buf.iv, &mut buf.pre);
        for r in 0..nch {
            let (x1, y1) = (buf.sx[r], buf.sy[r]);
            let l = m.mul(m.sub(nhy, y1), buf.iv[r]);
            let x3 = m.sub(m.sub(m.mul(l, l), x1), hx);
            buf.sx[r] = x3;
            buf.sy[r] = m.sub(m.mul(l, m.sub(x1, x3)), y1);
        }
    }
    out.retain(|&k| k >= -b2 && k <= b2);
    out.sort();
    out.dedup();
    Some(out)
}

/// All k in [-b2, b2] with (p+1-k) P = O, for the affine point P = (px, py) on y^2 = x^3 + a x + b.
fn kset(c: &Curve, p: u64, b2: i64, px: u64, py: u64, buf: &mut Buf) -> Vec<i64> {
    let w = (2 * b2 + 1) as u64;
    let mm = ((w as f64 / 2.0).sqrt().ceil() as u64).max(1); // baby steps 0..=mm
    let h = 2 * mm + 1;
    // baby steps: x(jP), j = 1..=mm
    buf.jac.clear();
    let mut q = Jac { x: px, y: py, z: c.m.one };
    for _ in 1..=mm {
        buf.jac.push(q);
        q = c.madd(q, px, py);
    }
    let (jac, aff, scr) = (&buf.jac, &mut buf.aff, &mut buf.scr);
    c.normalize(jac, aff, scr);
    buf.table_reset(mm as usize);
    let mut baby_inf = false; // some jP = O for j <= mm (small order)
    for j in 0..buf.aff.len() {
        match buf.aff[j] {
            // a repeated x means jP = ±j'P, i.e. order <= 2 mm: candidates could be lost, so reject the point
            Some((x, y)) => { if !buf.insert(x, y, j as u32 + 1) { baby_inf = true; } }
            None => baby_inf = true,
        }
    }
    if baby_inf { return vec![i64::MIN]; } // order <= mm: useless point, caller tries another
    // H = h P (affine), giant points G_i = (p + 1 - c0 - i h) P, c0 = -b2 + mm
    let hp = c.madd(c.dbl(q), px, c.m.neg(py)); // 2 (mm+1) P - P = (2 mm + 1) P = h P
    buf.jac.clear();
    buf.jac.push(hp);
    c.normalize(&buf.jac, &mut buf.aff, &mut buf.scr);
    let (hx, hy) = match buf.aff[0] { Some(v) => v, None => return vec![i64::MIN] };
    let nhy = c.m.neg(hy);
    let c0 = -b2 + mm as i64;
    let start = (p as i64 + 1 - c0) as u64;
    let mut g = c.scalar(start, px, py);
    let steps = ((2 * b2) as u64 / h + 2) as usize;
    let mut out = Vec::new();
    let batch = 64;
    let mut i0 = 0usize;
    while i0 < steps {
        let n = batch.min(steps - i0);
        buf.jac.clear();
        for _ in 0..n {
            buf.jac.push(g);
            g = c.madd(g, hx, nhy);
        }
        c.normalize(&buf.jac, &mut buf.aff, &mut buf.scr);
        for t in 0..n {
            let i = (i0 + t) as i64;
            let base = c0 + i * h as i64; // G_i = (p+1-base) P
            match buf.aff[t] {
                None => out.push(base),
                Some((x, y)) => {
                    if let Some((j, yj)) = buf.get(x) {
                        // G_i = (p+1-base) P equals  jP (same y) -> k = base + j,  or  -jP -> k = base - j
                        out.push(if y == yj { base + j as i64 } else { base - j as i64 });
                    }
                }
            }
        }
        i0 += n;
    }
    out.retain(|&k| k >= -b2 && k <= b2);
    out.sort();
    out.dedup();
    // verify each candidate exactly: (p+1-k) P == O (guards against x-collisions of distinct multiples)
    if cfg!(debug_assertions) { out.retain(|&k| c.scalar((p as i64 + 1 - k) as u64, px, py).z == 0); }
    out
}

fn ap_good(p: u64, a_big: &Big, b_big: &Big, buf: &mut Buf) -> i64 {
    let m = Mont::new(p);
    let a = m.to(a_big.modp(p));
    let b = m.to(b_big.modp(p));
    let b2 = isqrt(4 * p) as i64; // |a_p| <= 2 sqrt p  <=>  a_p^2 <= 4p
    // a_p = p + 1 - #E is even  <=>  #E even  <=>  E has a rational 2-torsion point  <=>  the cubic has a root
    let par: i64 = if cubic_has_root(&m, a, b) { 0 } else { 1 };
    let mut cand: Option<Vec<i64>> = None;
    let mut x0 = 1u64;
    loop {
        x0 += 1;
        if x0 > 10_000 { panic!("a_p not determined for p = {p}"); }
        let xm = m.to(x0);
        let d = m.add(m.mul(m.add(m.mul(xm, xm), a), xm), b); // x^3 + a x + b
        if d == 0 { continue; }
        let chi = jacobi(m.from(d), p) as i64;
        let d2 = m.mul(d, d);
        let curve = Curve { m: &m, a: m.mul(a, d2) };
        let (qx, qy) = (m.mul(xm, d), d2);
        let ks = match kset_fast(&curve, p, b2, qx, qy, par, 2, buf) {
            Some(v) => { buf.n_fast += 1; v }
            None => { buf.n_slow += 1; kset(&curve, p, b2, qx, qy, buf) }
        };
        if ks.first() == Some(&i64::MIN) || ks.is_empty() { continue; }
        let aps: Vec<i64> = ks.iter().map(|k| chi * k).collect();
        let next: Vec<i64> = match &cand { None => aps, Some(c) => c.iter().copied().filter(|v| aps.contains(v)).collect() };
        if next.is_empty() { panic!("inconsistent candidates for p = {p}"); }
        if next.len() == 1 { return next[0]; }
        cand = Some(next);
    }
}

/// Primes in [lo, hi) by segmented sieve.
fn primes_in(lo: u64, hi: u64) -> Vec<u64> {
    let lo = lo.max(2);
    if hi <= lo { return vec![]; }
    let r = isqrt(hi) + 1;
    let mut small = vec![true; (r + 1) as usize];
    let mut base = Vec::new();
    for i in 2..=r as usize {
        if small[i] { base.push(i as u64); let mut j = i * i; while j <= r as usize { small[j] = false; j += i; } }
    }
    let mut seg = vec![true; (hi - lo) as usize];
    for &q in &base {
        let mut s = ((lo + q - 1) / q * q).max(q * q);
        while s < hi { seg[(s - lo) as usize] = false; s += q; }
    }
    (lo..hi).filter(|&n| seg[(n - lo) as usize]).collect()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() == 5 && args[1] == "bench" {
        // ap bench <c4> <c6> <p>: per-component timings at one prime
        let p: u64 = primes_in(args[4].parse::<f64>().unwrap() as u64, args[4].parse::<f64>().unwrap() as u64 + 1000)[0];
        eprintln!("p = {p}");
        let a_big = Big::parse(&mul_small_str(&args[2], -27));
        let b_big = Big::parse(&mul_small_str(&args[3], -54));
        let m = Mont::new(p);
        let n = 20000u32;
        let tm = |name: &str, f: &mut dyn FnMut() -> u64| {
            let t = std::time::Instant::now(); let mut acc = 0u64;
            for _ in 0..n { acc = acc.wrapping_add(f()); }
            eprintln!("{name:12} {:8.1} ns  ({acc})", t.elapsed().as_nanos() as f64 / n as f64);
        };
        let mut x = 12345u64;
        tm("mont mul", &mut || { x = m.mul(x, x ^ 7); x });
        tm("inv", &mut || { x = m.inv(x | 1); x });
        tm("big_mod A+B", &mut || a_big.modp(p) + b_big.modp(p));
        tm("Mont::new", &mut || Mont::new(p).r2);
        let a = m.to(a_big.modp(p)); let b = m.to(b_big.modp(p));
        let xm = m.to(3); let d = m.add(m.mul(m.add(m.mul(xm, xm), a), xm), b);
        let d2 = m.mul(d, d); let cv = Curve { m: &m, a: m.mul(a, d2) };
        let (qx, qy) = (m.mul(xm, d), d2);
        tm("scalar p", &mut || cv.scalar(p + 1, qx, qy).x);
        let b2 = isqrt(4 * p) as i64;
        let mut buf = Buf::default();
        tm("cubic root?", &mut || cubic_has_root(&m, a, b) as u64);
        tm("kset_fast/2", &mut || kset_fast(&cv, p, b2, qx, qy, 0, 2, &mut buf).map(|v| v.len() as u64).unwrap_or(999));
        tm("kset_fast/1", &mut || kset_fast(&cv, p, b2, qx, qy, 0, 1, &mut buf).map(|v| v.len() as u64).unwrap_or(999));
        tm("kset", &mut || kset(&cv, p, b2, qx, qy, &mut buf).len() as u64);
        let ps = primes_in(p, p + 600_000); let mut it = 0usize; let mut buf2 = Buf::default();
        tm("ap_good", &mut || { it = (it + 1) % ps.len(); ap_good(ps[it], &a_big, &b_big, &mut buf2) as u64 });
        buf2.n_fast = 0; buf2.n_slow = 0;
        tm("ap_good", &mut || { it = (it + 1) % ps.len(); ap_good(ps[it], &a_big, &b_big, &mut buf2) as u64 });
        eprintln!("points fast {} slow {}", buf2.n_fast, buf2.n_slow);
        return;
    }
    if args.len() < 7 {
        eprintln!("usage: ap <c4> <c6> <small_ap.txt> <lo> <hi> <out_prefix> [threads]");
        std::process::exit(2);
    }
    let c4 = &args[1];
    let c6 = &args[2];
    // y^2 = x^3 - 27 c4 x - 54 c6  (isomorphic to E over F_p for p >= 5 good)
    let a_str = Arc::new(Big::parse(&mul_small_str(c4, -27)));
    let b_str = Arc::new(Big::parse(&mul_small_str(c6, -54)));
    let mut table: HashMap<u64, i64> = HashMap::new();
    for line in std::fs::read_to_string(&args[3]).unwrap().lines() {
        let mut it = line.split_whitespace();
        if let (Some(p), Some(a)) = (it.next(), it.next()) { table.insert(p.parse().unwrap(), a.parse().unwrap()); }
    }
    let table = Arc::new(table);
    let lo: u64 = args[4].parse::<f64>().unwrap() as u64;
    let hi: u64 = args[5].parse::<f64>().unwrap() as u64;
    let prefix = args[6].clone();
    let threads: usize = args.get(7).map(|s| s.parse().unwrap())
        .unwrap_or_else(|| std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1));
    let blk: u64 = std::env::var("AP_BLOCK").ok().map(|v| v.parse::<f64>().unwrap() as u64).unwrap_or(100_000_000);
    let legacy = std::env::var("AP_LEGACY").is_ok(); // also write .p.u64/.ap.i64 (for comparisons)
    let nblk = ((hi - lo) + blk - 1) / blk;
    let next = Arc::new(AtomicUsize::new(0));
    let t0 = std::time::Instant::now();
    let handles: Vec<_> = (0..threads).map(|_| {
        let (next, table, a_str, b_str, prefix) = (next.clone(), table.clone(), a_str.clone(), b_str.clone(), prefix.clone());
        let legacy = legacy;
        std::thread::spawn(move || {
            let mut buf = Buf::default();
            loop {
                let i = next.fetch_add(1, Ordering::SeqCst) as u64;
                if i >= nblk {
                    FAST.fetch_add(buf.n_fast, Ordering::Relaxed);
                    SLOW.fetch_add(buf.n_slow, Ordering::Relaxed);
                    break;
                }
                let (b_lo, b_hi) = (lo + i * blk, (lo + (i + 1) * blk).min(hi));
                let name = format!("{prefix}.{b_lo:012}-{b_hi:012}");
                if std::path::Path::new(&format!("{name}.ok")).exists() { continue; }
                let ps = primes_in(b_lo, b_hi);
                let mut aps = Vec::with_capacity(ps.len());
                for &p in &ps {
                    let ap = if p < P_SMALL || table.contains_key(&p) {
                        *table.get(&p).unwrap_or_else(|| panic!("small/bad prime {p} missing from table"))
                    } else {
                        ap_good(p, &a_str, &b_str, &mut buf)
                    };
                    aps.push(ap);
                }
                write_rap(&format!("{name}.rap"), b_lo, b_hi, &ps, &aps);
                if legacy {
                    let mut fp = std::fs::File::create(format!("{name}.p.u64")).unwrap();
                    for p in &ps { fp.write_all(&p.to_le_bytes()).unwrap(); }
                    let mut fa = std::fs::File::create(format!("{name}.ap.i64")).unwrap();
                    for a in &aps { fa.write_all(&a.to_le_bytes()).unwrap(); }
                }
                std::fs::write(format!("{name}.ok"), format!("{}\n", ps.len())).unwrap();
                eprintln!("{name}: {} primes, {:.1?} elapsed", ps.len(), t0.elapsed());
            }
        })
    }).collect();
    for h in handles { h.join().unwrap(); }
    eprintln!("points: fast {} / fallback {}", FAST.load(Ordering::Relaxed), SLOW.load(Ordering::Relaxed));
}

/// decimal string times small signed integer, as decimal string
fn mul_small_str(s: &str, k: i64) -> String {
    let (neg, digits) = match s.strip_prefix('-') { Some(d) => (true, d), None => (false, s) };
    let neg = neg ^ (k < 0);
    let k = k.unsigned_abs() as u128;
    let mut out: Vec<u8> = Vec::new();
    let mut carry: u128 = 0;
    for &ch in digits.as_bytes().iter().rev() {
        let v = (ch - b'0') as u128 * k + carry;
        out.push(b'0' + (v % 10) as u8);
        carry = v / 10;
    }
    while carry > 0 { out.push(b'0' + (carry % 10) as u8); carry /= 10; }
    while out.len() > 1 && *out.last().unwrap() == b'0' { out.pop(); }
    out.reverse();
    let body = String::from_utf8(out).unwrap();
    if neg && body != "0" { format!("-{body}") } else { body }
}


/// Compact block file: b"RAP1", lo u64, hi u64, n u64, gaps u16[n] (p_0 - lo, p_i - p_{i-1}), a_p i32[n]; little-endian.
fn write_rap(path: &str, lo: u64, hi: u64, ps: &[u64], aps: &[i64]) {
    let mut v: Vec<u8> = Vec::with_capacity(32 + 6 * ps.len());
    v.extend_from_slice(b"RAP1");
    for x in [lo, hi, ps.len() as u64] { v.extend_from_slice(&x.to_le_bytes()); }
    let mut prev = lo;
    for &p in ps {
        let g = p - prev;
        assert!(g < 65536, "prime gap too large for u16");
        v.extend_from_slice(&(g as u16).to_le_bytes());
        prev = p;
    }
    for &a in aps { v.extend_from_slice(&i32::try_from(a).unwrap().to_le_bytes()); }
    let tmp = format!("{path}.tmp");
    std::fs::write(&tmp, &v).unwrap();
    std::fs::rename(&tmp, path).unwrap();
}
