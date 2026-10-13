import json, math, sys, numpy as np
sys.path.insert(0,'/home/claude/highrank/radix-rank/scripts'); import certify as C
from flint import arb, fmpq, fmpq_poly
def cert_mem(n, key):
    c=json.load(open(n+'_ap.json')); d=json.loads(open(n+'_scan.jsonl').readline()); h=d[key]
    cv=dict(N=int(c['N_str']),primes=c['primes'],ap=c['ap'])
    b=C.certify(cv,repr(h['l']),h['coeffs']); t=d['rank']+(2 if key.endswith('r2') else 1)
    return h, b, t
def cert_bin(binjson, scanfile, key, fixedL=None):
    c=json.loads(open(binjson).readline()); d=json.loads(open(scanfile).readline())
    h=d[key] if key else dict(l=fixedL,primes=None,coeffs=[1/(2*fixedL)]+[0.0]*9)
    L=h['l']; X=math.exp(2*L)*1.000001
    p=np.fromfile(c['bin_p'],dtype='<u8'); a=np.fromfile(c['bin_ap'],dtype='<i8'); k=int(np.searchsorted(p,X))+1
    bad=set(c['bad'])
    # streamed prime sum (same maths as certify.py)
    g=C.g_poly(h['coeffs'],repr(L)); gc=[fmpq(x) for x in g.coeffs()]; g0=arb(g(fmpq(0)))
    W=g0*(arb(int(c['N_str'])).log()-2*(2*arb.pi()).log())+C.arch_term(g,repr(L)); twoL=2*float(C.fq(repr(L))); S=arb(0)
    for pp,ap in zip(p[:k].tolist(),a[:k].tolist()):
        lp=arb(pp).log()
        if float(lp.lower())>=twoL: break
        s1=arb(ap)/arb(pp).sqrt(); sp,sc,kk=arb(2),s1,1
        while True:
            ln=kk*lp
            if float(ln.lower())>=twoL: break
            if float(ln.upper())>=twoL: raise ValueError('edge')
            ck=s1**kk if pp in bad else sc
            S+=lp*ck/arb(pp)**(arb(kk)/2)*C.horner_arb(gc,ln)
            if pp not in bad: sp,sc=sc,s1*sc-sp
            kk+=1
    W-=2*S
    f=sum((C.fq(x)*fmpq_poly.legendre_p(2*i) for i,x in enumerate(h['coeffs'])),fmpq_poly(0)); F=f.integral()
    return h, W/arb(C.fq(repr(L))*(F(fmpq(1))-F(fmpq(-1))))**2
for n in ('E20','E21','E22','E23','E24'):
    for key in ('std_r2','radix_r2','radix_r1'):
        d=json.loads(open(n+'_scan.jsonl').readline())
        if not d[key]: continue
        h,b,t=cert_mem(n,key)
        print(f"{n} {key:9s} L={h['l']:<5} primes={h['primes']:>8}  Arb {b.str(10,radius=True)} < {t}: {'PROVEN' if b.upper()<t else 'NOT'}",flush=True)
for key in ('std_r2','radix_r2'):
    h,b=cert_bin('E27_scan.json','E27_scan.jsonl',key)
    print(f"E27 {key:9s} L={h['l']:<5} primes={h['primes']:>8}  Arb {b.str(10,radius=True)} < 29: {'PROVEN' if b.upper()<29 else 'NOT'}",flush=True)
