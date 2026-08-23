#!/usr/bin/env python3
"""Exact-rational interval certificate for the near-aligned JNR theorem.

The polynomial kernel and the initial Bernstein coefficients are constructed
exactly over Q.  De Casteljau subdivision is evaluated in binary64 only as a
center representation.  Every center is interpreted exactly as a dyadic
rational and is accompanied by an *exact Fraction error radius* that provably
encloses the corresponding exact Bernstein coefficient.  Error radii are
propagated by a rational recurrence.  Every pruning decision is therefore an
exact rational interval sign comparison; no floating-point inequality is used
to establish emptiness.

This is a rigorous interval-arithmetic certificate with exact rational
endpoints for every sign decision.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from math import comb
import argparse
import hashlib
import json
import math
import time

import numpy as np
import sympy as sp

p, a, b, y, u, v = sp.symbols('p a b y u v', real=True)
VARS = (p, a, b, y, u, v)
VAR_NAMES = ('p','a','b','y','u','v')
M = 1-p**2
Sr = (1+p)*a
Si = (1-p)*b
qr = u**2-v**2-Sr*u+Si*v+p
qi = 2*u*v-Sr*v-Si*u
G = sp.expand(qr**2+qi**2 + y*(y + sp.Rational(1,2)*((2*u-Sr)**2+(2*v-Si)**2)
                                      +1+p**2-sp.Rational(1,2)*(Sr**2+Si**2)))
D0 = 1+y+u**2+v**2-2*(a*u+b*v)
A0 = 1-y-u**2-v**2
L2 = sp.expand(A0**2+4*y)
H = sp.expand(A0*G-2*M*y*D0)
Gy,Gu,Gv = sp.diff(G,y),sp.diff(G,u),sp.diff(G,v)
Hy,Hu,Hv = sp.diff(H,y),sp.diff(H,u),sp.diff(H,v)
Cu = sp.expand(Gy*Hu-Gu*Hy)
Cv = sp.expand(Gy*Hv-Gv*Hy)
Ds = sp.expand(M*(1-a**2-b**2))
Qs = sp.expand(3+p**2-(1+p)**2*a**2-(1-p)**2*b**2)
Ktilde = sp.expand(M*Ds*H-Qs*G**2)
DOMAIN = (Fraction(0),Fraction(1,2), Fraction(-1),Fraction(1), Fraction(-1),Fraction(1),
          Fraction(0),Fraction(1), Fraction(-1),Fraction(1), Fraction(-1),Fraction(1))

# Four times binary64 unit roundoff.  This is an exact rational constant.
UROUND = Fraction(1, 2**51)


def F(q):
    q = sp.Rational(q)
    return Fraction(int(q.p), int(q.q))


def exact_bernstein_dict(expr):
    """Tensor-product Bernstein coefficients on DOMAIN, exactly in Q."""
    P = sp.Poly(expr,*VARS)
    deg = [P.degree(x) for x in VARS]
    D = {mon:F(co) for mon,co in P.terms()}
    for ax,(lo,hi,n) in enumerate(zip(DOMAIN[::2],DOMAIN[1::2],deg)):
        w = hi-lo
        T = [[Fraction(0) for _ in range(n+1)] for __ in range(n+1)]
        for k in range(n+1):
            for m in range(n+1):
                sm = Fraction(0)
                for i in range(min(m,k)+1):
                    sm += Fraction(comb(m,i)*comb(k,i),comb(n,i))*lo**(m-i)*w**i
                T[k][m] = sm
        groups = defaultdict(lambda:[Fraction(0) for _ in range(n+1)])
        for mon,c in D.items():
            key = mon[:ax]+mon[ax+1:]
            groups[key][mon[ax]] = c
        ND = {}
        for key,vec in groups.items():
            for k in range(n+1):
                z = sum((T[k][m]*vec[m] for m in range(n+1)),Fraction(0))
                if z:
                    ND[key[:ax]+(k,)+key[ax:]] = z
        D = ND
    return D,deg


def next_power_two_strict(q: Fraction) -> int:
    """Integer power of two strictly larger than q>=0."""
    if q <= 0:
        return 1
    n = 1
    # q sizes here are modest; exact loop keeps this dependency-free.
    while Fraction(n) <= q:
        n <<= 1
    return n


def exact_to_float_tensor(expr):
    """Exact initial coefficients -> float centers + exact rational radius."""
    t0 = time.time()
    D,deg = exact_bernstein_dict(expr)
    shape = tuple(d+1 for d in deg)
    B = np.zeros(shape,dtype=np.float64)
    maxerr = Fraction(0)
    exact_maxabs = Fraction(0)
    for mon,q in D.items():
        f = float(q)
        if not math.isfinite(f):
            raise ArithmeticError('non-finite initial Bernstein center')
        B[mon] = f
        err = abs(q-Fraction.from_float(f))
        maxerr = max(maxerr,err)
        exact_maxabs = max(exact_maxabs,abs(q))
    # True coefficients remain within the initial exact convex hull under
    # de Casteljau.  Choose a generous exact power-of-two bound on all stored
    # centers.  The propagated radii stay <1 in these runs, verified at every
    # split, so 2*(max_exact+1) leaves a strict safety margin.
    center_bound = next_power_two_strict(2*(exact_maxabs+1))
    # Close the root-box center-bound invariant explicitly.  np.max only
    # selects a stored binary64 value; Fraction.from_float converts that
    # selected center to its exact dyadic rational before comparison.
    root_max = float(np.max(np.abs(B)))
    if Fraction.from_float(root_max) >= Fraction(center_bound):
        raise ArithmeticError('initial center magnitude exceeded rigorous bound')
    if maxerr >= 1:
        raise ArithmeticError('initial rational error radius reached 1')
    return B,deg,maxerr,Fraction(center_bound),{
        'seconds':time.time()-t0,
        'degree':deg,
        'shape':list(shape),
        'nonzero_exact_bernstein_coefficients':len(D),
        'exact_initial_max_abs':str(exact_maxabs),
        'exact_initial_to_binary64_radius':str(maxerr),
        'center_abs_bound':str(center_bound),
    }


def split_float_bern(B,axis):
    X = np.moveaxis(B,axis,0)
    n = X.shape[0]-1
    temp = X.copy(); L = np.empty_like(X); R = np.empty_like(X)
    L[0] = temp[0]; R[n] = temp[n]
    for r in range(1,n+1):
        temp = (temp[:-1]+temp[1:])*0.5
        L[r] = temp[0]; R[n-r] = temp[-1]
    if not np.all(np.isfinite(L)) or not np.all(np.isfinite(R)):
        raise ArithmeticError('non-finite de Casteljau center')
    return np.moveaxis(L,0,axis),np.moveaxis(R,0,axis)


def expr_hash(expr):
    return hashlib.sha256(str(sp.Poly(expr,*VARS).as_dict()).encode()).hexdigest()


class PolyState:
    __slots__=('name','B','sense','deg','err','center_bound')
    def __init__(self,name,B,sense,deg,err,center_bound):
        self.name=name; self.B=B; self.sense=sense; self.deg=tuple(deg)
        self.err=Fraction(err); self.center_bound=Fraction(center_bound)

    def split(self,axis):
        L,R = split_float_bern(self.B,axis)
        # A child coefficient is formed by at most deg[axis] rounded averages.
        # For one average fl((x+y)/2), if both input centers have magnitude <B,
        # the additional absolute rounding error is safely <= 2^-51 B.
        # Existing endpoint uncertainty averages without amplification.
        e = self.err + self.deg[axis]*UROUND*self.center_bound
        if e >= 1:
            raise ArithmeticError('error radius reached 1; center-bound invariant not closed')
        # Stored centers are bounded by exact coefficient max + e < bound by construction.
        maxLR = max(float(np.max(np.abs(L))),float(np.max(np.abs(R))))
        if Fraction.from_float(maxLR) >= self.center_bound:
            raise ArithmeticError('center magnitude exceeded rigorous bound')
        return (PolyState(self.name,L,self.sense,self.deg,e,self.center_bound),
                PolyState(self.name,R,self.sense,self.deg,e,self.center_bound))

    def prune(self):
        # np.min/max only select an existing stored binary64 center.  Convert it
        # exactly to its dyadic Fraction and compare the rational interval
        # endpoint with zero.  The sign decision itself is exact arithmetic.
        lohat = Fraction.from_float(float(np.min(self.B)))
        hihat = Fraction.from_float(float(np.max(self.B)))
        if self.sense=='zero':
            lower = lohat-self.err
            upper = hihat+self.err
            if lower>0: return True,lower,'positive'
            if upper<0: return True,-upper,'negative'
        elif self.sense=='le0':
            lower = lohat-self.err
            if lower>0: return True,lower,'positive'
        elif self.sense=='ge0':
            upper = hihat+self.err
            if upper<0: return True,-upper,'negative'
        return False,None,None


def build_poly(name,expr,sense):
    B,deg,e0,bound,info = exact_to_float_tensor(expr)
    return PolyState(name,B,sense,deg,e0,bound),info


def copy_state(st):
    return PolyState(st.name,st.B.copy(),st.sense,st.deg,st.err,st.center_bound)


def run_certificate(delta,h0,out_json,common=None,max_boxes=250000):
    delta=Fraction(delta); h0=Fraction(h0)
    X=(Fraction(8)+3*delta)/(Fraction(2)+3*delta)
    Cap=sp.expand(Qs-sp.Rational(X.numerator,X.denominator)*Ds)
    Nh=sp.expand(H*Qs-sp.Rational(h0.numerator,h0.denominator)**2*M*Ds*L2)
    build={}; states=[]
    if common is None:
        for name,E,sense in [('Ktilde',Ktilde,'zero'),('Cu',Cu,'zero'),('Cv',Cv,'zero')]:
            st,info=build_poly(name,E,sense); states.append(st); build[name]=info
    else:
        for st,info in common:
            states.append(copy_state(st)); build[st.name]=dict(info)
    for name,E,sense in [('Nh',Nh,'ge0'),('Cap',Cap,'le0')]:
        st,info=build_poly(name,E,sense); states.append(st); build[name]=info

    order=[4,3,0,1,2] # Cap, Nh, Ktilde, Cu, Cv
    stack=[((0,0,0,0,0,0),states)]
    visited=0; prunes=Counter()
    minmargin={s.name:None for s in states}
    maxerr_seen={s.name:s.err for s in states}
    max_axis_depth=[0]*6; deepest=(0,0,0,0,0,0); t0=time.time()
    while stack:
        depths,S=stack.pop(); visited+=1
        if sum(depths)>sum(deepest): deepest=depths
        for i,d in enumerate(depths): max_axis_depth[i]=max(max_axis_depth[i],d)
        killed=False
        for idx in order:
            st=S[idx]
            yes,margin,_=st.prune()
            maxerr_seen[st.name]=max(maxerr_seen[st.name],st.err)
            if yes:
                prunes[st.name]+=1
                if minmargin[st.name] is None or margin<minmargin[st.name]: minmargin[st.name]=margin
                killed=True; break
        if killed: continue
        axis=min(range(6),key=lambda i:depths[i])
        nd=list(depths); nd[axis]+=1; nd=tuple(nd)
        LS=[]; RS=[]
        for st in S:
            l,r=st.split(axis); LS.append(l); RS.append(r)
        stack.append((nd,RS)); stack.append((nd,LS))
        if visited>max_boxes: raise RuntimeError('did not close')
    elapsed=time.time()-t0
    eps_exact=(Fraction(2,3)+delta)*(1-h0)/(1+h0)
    rho2=1-h0*h0
    def ratrec(q):
        return None if q is None else {'exact':str(q),'decimal':float(q)}
    result={
      'status':'CERTIFIED_EMPTY',
      'delta':str(delta),'X_delta':str(X),'h0':str(h0),
      'rho_delta_squared_lower_bound':str(rho2),'epsilon_delta_lower_bound':str(eps_exact),
      'outer_domain':{'p':['0','1/2'],'a':['-1','1'],'b':['-1','1'],'y':['0','1'],'u':['-1','1'],'v':['-1','1']},
      'polynomial_system':{
          'Ktilde':'M Ds H - Qs G^2 = 0',
          'Cu':'G_y H_u - G_u H_y = 0',
          'Cv':'G_y H_v - G_v H_y = 0',
          'Cap':'Qs - X_delta Ds <= 0',
          'Nh':'H Qs - h0^2 M Ds L^2 >= 0'},
      'statement':'No simultaneous solution of Ktilde=Cu=Cv=0, Cap<=0 and Nh>=0 exists on the outer box. Therefore no physical ordinary node with Xi_J<=X_delta and h>=h0 exists.',
      'visited_boxes':visited,'prune_counts':dict(prunes),
      'max_axis_depths':dict(zip(VAR_NAMES,max_axis_depth)),'deepest_depth_vector':list(deepest),
      'min_certified_rational_interval_sign_margin':{k:ratrec(q) for k,q in minmargin.items()},
      'max_exact_rational_error_radius_seen':{k:ratrec(q) for k,q in maxerr_seen.items()},
      'build':build,'subdivision_seconds':elapsed,
      'arithmetic_note':('Initial tensor-product Bernstein coefficients are exact fractions. Binary64 values are only centers. '
         'Each center is interpreted exactly as a dyadic Fraction, with an exact rational radius. Under a half-split along '
         'axis i the rigorous uniform radius is E_child = E_parent + d_i*2^-51*B, with B an exact power-of-two bound on '
         'all stored centers; the bound invariant is checked at every split. Every pruning sign is decided by an exact '
         'Fraction comparison of the interval endpoint with zero. No floating-point inequality or approximate error-envelope '
         'comparison enters the proof.'),
      'bernstein_note':'The exact polynomial range on a box lies in the convex hull of its exact Bernstein coefficients; each such coefficient lies in the exact rational interval centered on the stored dyadic value.',
      'hashes':{'Ktilde':expr_hash(Ktilde),'Cu':expr_hash(Cu),'Cv':expr_hash(Cv),'Nh':expr_hash(Nh),'Cap':expr_hash(Cap)}
    }
    with open(out_json,'w') as f: json.dump(result,f,indent=2,sort_keys=True)
    return result


def build_common():
    out=[]
    for name,E,sense in [('Ktilde',Ktilde,'zero'),('Cu',Cu,'zero'),('Cv',Cv,'zero')]:
        st,info=build_poly(name,E,sense); out.append((st,info)); print('built',name,info,flush=True)
    return out


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--case',choices=['0.02','0.05','both'],default='both'); ap.add_argument('--prefix',default='jnr_rational_interval_v3'); args=ap.parse_args()
    common=build_common(); cases=[]
    if args.case in ('0.02','both'): cases.append((Fraction(1,50),Fraction(24,25),f'{args.prefix}_delta002.json'))
    if args.case in ('0.05','both'): cases.append((Fraction(1,20),Fraction(9,10),f'{args.prefix}_delta005.json'))
    for delta,h0,fn in cases:
        print('RUN',delta,h0,flush=True); R=run_certificate(delta,h0,fn,common)
        keys=['status','delta','X_delta','h0','rho_delta_squared_lower_bound','epsilon_delta_lower_bound','visited_boxes','prune_counts','max_axis_depths','min_certified_rational_interval_sign_margin','max_exact_rational_error_radius_seen','subdivision_seconds']
        print(json.dumps({k:R[k] for k in keys},indent=2),flush=True)

if __name__=='__main__': main()
