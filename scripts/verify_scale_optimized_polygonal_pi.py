#!/usr/bin/env python3
"""Verify scale-b polygonal acceleration, stability, and constructibility checks."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import mpmath as mp

BRANCHES = (2,3,6,9,27)
ORDERS = (4,8,12,24)
SLOPE_N = (16,20,24,32)

def polygon(n: int) -> mp.mpf:
    x = mp.mpf(n)
    return x*mp.sin(mp.pi/x)

def weights(b: int, r: int) -> list[mp.mpf]:
    xs=[mp.mpf(b)**(-2*m) for m in range(r)]
    out=[]
    for m,xm in enumerate(xs):
        w=mp.mpf(1)
        for i,xi in enumerate(xs):
            if i != m:
                w *= (-xi)/(xm-xi)
        out.append(w)
    return out

def accelerated(N: int,b: int,r: int) -> mp.mpf:
    w=weights(b,r)
    return mp.fsum(w[m]*polygon((b**m)*N) for m in range(r))

def stability_constant(b: int) -> mp.mpf:
    q=mp.mpf(b)**-2
    prod=mp.mpf(1)
    ell=1
    threshold=mp.power(10, -(mp.mp.dps-30))
    while q**ell > threshold:
        t=q**ell
        prod *= (1+t)/(1-t)
        ell += 1
    return prod

def regression_slope(xs,ys):
    xb=mp.fsum(xs)/len(xs); yb=mp.fsum(ys)/len(ys)
    return mp.fsum((x-xb)*(y-yb) for x,y in zip(xs,ys))/mp.fsum((x-xb)**2 for x in xs)

def observed_slope(b: int,r: int) -> mp.mpf:
    xs=[mp.log10(N) for N in SLOPE_N]
    ys=[mp.log10(abs(accelerated(N,b,r)-mp.pi)) for N in SLOPE_N]
    return regression_slope(xs,ys)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--precision',type=int,default=1300)
    p.add_argument('--output-dir',type=Path,default=Path('output'))
    a=p.parse_args()
    if a.precision < 1100:
        raise ValueError('Use at least 1100 decimal digits for the r=24, b=27 audit.')
    mp.mp.dps=a.precision
    a.output_dir.mkdir(parents=True,exist_ok=True)

    rows=[]; moments_ok=True; slopes_ok=True; bounds_ok=True
    for b in BRANCHES:
        Sb=stability_constant(b)
        for r in ORDERS:
            w=weights(b,r)
            unity=abs(mp.fsum(w)-1)
            moment=max([abs(mp.fsum(w[m]*mp.mpf(b)**(-2*j*m) for m in range(r))) for j in range(1,r)] or [mp.mpf(0)])
            moments_ok &= unity < mp.mpf('1e-900') and moment < mp.mpf('1e-900')
            slope=observed_slope(b,r)
            slopes_ok &= abs(slope+2*r) < mp.mpf('5e-3')
            N=16
            err=abs(accelerated(N,b,r)-mp.pi)
            bound=Sb*mp.pi**(2*r+1)/mp.factorial(2*r+1)*mp.mpf(N)**(-2*r)
            bounds_ok &= err <= bound
            rows.append({
                'b':b,'r':r,'stability_constant':mp.nstr(Sb,35),
                'weight_l1':mp.nstr(mp.fsum(abs(x) for x in w),35),
                'unity_residual':mp.nstr(unity,8),'max_moment_residual':mp.nstr(moment,8),
                'observed_slope':mp.nstr(slope,30),'predicted_slope':-2*r,
                'N_for_bound':N,'actual_error':mp.nstr(err,30),'analytic_bound':mp.nstr(bound,30)
            })
    with (a.output_dir/'scale_acceleration_audit.csv').open('w',newline='',encoding='utf-8') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)

    # dyadic chord identity
    chord_ok=True
    c=mp.mpf(1); N=6
    for _ in range(8):
        cnext=c/mp.sqrt(2+mp.sqrt(4-c*c))
        exact=2*mp.sin(mp.pi/(2*N))
        chord_ok &= abs(cnext-exact) < mp.mpf('1e-1100')
        c,N=cnext,2*N

    # linear tripling seed asymptotic
    target=-8*mp.pi**3/81
    vals=[]
    for N in (100,200,400,800,1600):
        cN=2*mp.sin(mp.pi/N); c3=2*mp.sin(mp.pi/(3*N))
        scaled=mp.mpf(N)**3*(cN/3-c3)
        vals.append({'N':N,'scaled_error':mp.nstr(scaled,35),'target':mp.nstr(target,35)})
    seed_ok=abs(mp.mpf(vals[-1]['scaled_error'])-target) < mp.mpf('1e-5')
    with (a.output_dir/'linear_tripling_seed_audit.csv').open('w',newline='',encoding='utf-8') as f:
        wr=csv.DictWriter(f,fieldnames=list(vals[0].keys())); wr.writeheader(); wr.writerows(vals)

    expected={2:'1.969260353668269',3:'1.285210586391831',6:'1.058822194456829',9:'1.025316406852818',27:'1.002751031629743'}
    constants_ok=True
    for b,t in expected.items():
        constants_ok &= abs(stability_constant(b)-mp.mpf(t)) < mp.mpf('1e-15')

    overall=bool(moments_ok and slopes_ok and bounds_ok and chord_ok and seed_ok and constants_ok)
    summary={
        'precision_dps':a.precision,
        'branches':list(BRANCHES),
        'orders':list(ORDERS),
        'moment_conditions_passed':bool(moments_ok),
        'stability_constants_passed':bool(constants_ok),
        'representative_slopes_passed':bool(slopes_ok),
        'analytic_remainder_bounds_passed':bool(bounds_ok),
        'dyadic_chord_recursion_passed':bool(chord_ok),
        'linear_tripling_seed_asymptotic_passed':bool(seed_ok),
        'verification_passed':overall,
    }
    (a.output_dir/'verification_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    (a.output_dir/'verification_summary.txt').write_text('\n'.join(f'{k}: {v}' for k,v in summary.items())+'\n',encoding='utf-8')
    print('\n'.join(f'{k}: {v}' for k,v in summary.items()))
    if not overall: raise SystemExit('Verification failed.')
if __name__=='__main__': main()
