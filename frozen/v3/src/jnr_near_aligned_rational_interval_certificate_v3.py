#!/usr/bin/env python3
from __future__ import annotations

from fractions import Fraction
from collections import defaultdict, Counter
from math import comb
import hashlib, json, time
import numpy as np
import sympy as sp

p,a,b,y,u,v = sp.symbols('p a b y u v', real=True)
VARS=(p,a,b,y,u,v)
M=1-p**2
Sr=(1+p)*a
Si=(1-p)*b
qr=u**2-v**2-Sr*u+Si*v+p
qi=2*u*v-Sr*v-Si*u
G=sp.expand(qr**2+qi**2 + y*(y + sp.Rational(1,2)*((2*u-Sr)**2+(2*v-Si)**2)
    +1+p**2-sp.Rational(1,2)*(Sr**2+Si**2)))
D0=1+y+u**2+v**2-2*(a*u+b*v)
A0=1-y-u**2-v**2
L2=sp.expand(A0**2+4*y)
H=sp.expand(A0*G-2*M*y*D0)
Gy,Gu,Gv=sp.diff(G,y),sp.diff(G,u),sp.diff(G,v)
Hy,Hu,Hv=sp.diff(H,y),sp.diff(H,u),sp.diff(H,v)
Cu=sp.expand(Gy*Hu-Gu*Hy)
Cv=sp.expand(Gy*Hv-Gv*Hy)
Ds=sp.expand(M*(1-a**2-b**2))
Qs=sp.expand(3+p**2-(1+p)**2*a**2-(1-p)**2*b**2)
Ktilde=sp.expand(M*Ds*H-Qs*G**2)
DOMAIN=(Fraction(0),Fraction(1,2), Fraction(-1),Fraction(1), Fraction(-1),Fraction(1),
        Fraction(0),Fraction(1), Fraction(-1),Fraction(1), Fraction(-1),Fraction(1))
UROUND_Q = Fraction(1, 2**51) # deliberately 4x binary64 unit roundoff


def F(q):
    q=sp.Rational(q)
    return Fraction(int(q.p),int(q.q))


def exact_bernstein_dict(expr):
    P=sp.Poly(expr,*VARS)
    deg=[P.degree(x) for x in VARS]
    D={mon:F(co) for mon,co in P.terms()}
    for ax,(lo,hi,n) in enumerate(zip(DOMAIN[::2],DOMAIN[1::2],deg)):
        w=hi-lo
        T=[[Fraction(0) for _ in range(n+1)] for __ in range(n+1)]
        for k in range(n+1):
            for m in range(n+1):
                sm=Fraction(0)
                for i in range(min(m,k)+1):
                    sm += Fraction(comb(m,i)*comb(k,i),comb(n,i))*lo**(m-i)*w**i
                T[k][m]=sm
        groups=defaultdict(lambda:[Fraction(0) for _ in range(n+1)])
        for mon,c in D.items():
            key=mon[:ax]+mon[ax+1:]
            groups[key][mon[ax]]=c
        ND={}
        for key,vec in groups.items():
            for k in range(n+1):
                z=sum((T[k][m]*vec[m] for m in range(n+1)),Fraction(0))
                if z:
                    mon=key[:ax]+(k,)+key[ax:]
                    ND[mon]=z
        D=ND
    return D,deg


def exact_to_float_tensor(expr):
    D,deg=exact_bernstein_dict(expr)
    shape=tuple(d+1 for d in deg)
    B=np.zeros(shape,dtype=np.float64)
    maxerr=Fraction(0)
    for mon,q in D.items():
        f=float(q)
        B[mon]=f
        err=abs(q-Fraction.from_float(f))
        if err>maxerr:
            maxerr=err
    maxabs=float(np.max(np.abs(B)))
    # Keep the same deliberately loose floating magnitude bound as the public engine,
    # but interpret the resulting dyadic exactly in all proof comparisons.
    Mbound_fp=np.nextafter(2.0*(maxabs+1.0),np.inf)
    Mbound_q=Fraction.from_float(float(Mbound_fp))
    return B,deg,maxerr,Mbound_q


def split_float_bern(B,axis):
    X=np.moveaxis(B,axis,0)
    n=X.shape[0]-1
    temp=X.copy(); L=np.empty_like(X); R=np.empty_like(X)
    L[0]=temp[0]; R[n]=temp[n]
    for r in range(1,n+1):
        temp=(temp[:-1]+temp[1:])*0.5
        L[r]=temp[0]; R[n-r]=temp[-1]
    return np.moveaxis(L,0,axis),np.moveaxis(R,0,axis)


def expr_hash(expr):
    return hashlib.sha256(str(sp.Poly(expr,*VARS).as_dict()).encode()).hexdigest()


class PolyState:
    __slots__=('name','B','sense','deg','err','Mbound')
    def __init__(self,name,B,sense,deg,err,Mbound):
        self.name=name; self.B=B; self.sense=sense; self.deg=tuple(deg); self.err=err; self.Mbound=Mbound
    def split(self,axis):
        L,R=split_float_bern(self.B,axis)
        # Exact rational error radius. At most deg[axis] rounded averages
        # contribute to any child coefficient.
        e=self.err + self.deg[axis]*UROUND_Q*self.Mbound
        return (PolyState(self.name,L,self.sense,self.deg,e,self.Mbound),
                PolyState(self.name,R,self.sense,self.deg,e,self.Mbound))
    def prune(self):
        # Floating comparisons are only a cheap prefilter.  A box is pruned
        # only after the selected dyadic center has been converted back to an
        # exact Fraction and compared with the exact rational error radius.
        lo=float(np.min(self.B)); hi=float(np.max(self.B)); ef=float(self.err)
        e=self.err
        if self.sense=='zero':
            if lo>ef:
                lo_q=Fraction.from_float(lo)
                if lo_q-e>0:
                    return True,lo_q-e,'positive'
            if hi<-ef:
                hi_q=Fraction.from_float(hi)
                if hi_q+e<0:
                    return True,-hi_q-e,'negative'
        elif self.sense=='le0':
            if lo>ef:
                lo_q=Fraction.from_float(lo)
                if lo_q-e>0:
                    return True,lo_q-e,'positive'
        elif self.sense=='ge0':
            if hi<-ef:
                hi_q=Fraction.from_float(hi)
                if hi_q+e<0:
                    return True,-hi_q-e,'negative'
        return False,None,None


def build_poly(name,expr,sense):
    t=time.time(); B,deg,e0,Mbound=exact_to_float_tensor(expr)
    return PolyState(name,B,sense,deg,e0,Mbound), {
        'seconds':time.time()-t,'degree':deg,'shape':list(B.shape),
        'max_initial_exact_to_float_error':str(e0),
        'max_initial_exact_to_float_error_float':float(e0),
        'roundoff_bound_initial':str(e0),
        'roundoff_bound_initial_float':float(e0),
        'Mbound':str(Mbound),'Mbound_float':float(Mbound)}


def run_certificate(delta,h0,out_json,common=None):
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
            states.append(PolyState(st.name,st.B.copy(),st.sense,st.deg,st.err,st.Mbound)); build[st.name]=dict(info)
    for name,E,sense in [('Nh',Nh,'ge0'),('Cap',Cap,'le0')]:
        st,info=build_poly(name,E,sense); states.append(st); build[name]=info
    # test Cap, Nh, Ktilde, Cu, Cv in that order
    order=[4,3,0,1,2]
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
            yes,margin,sgn=S[idx].prune()
            if S[idx].err>maxerr_seen[S[idx].name]: maxerr_seen[S[idx].name]=S[idx].err
            if yes:
                prunes[S[idx].name]+=1
                old=minmargin[S[idx].name]
                if old is None or margin<old: minmargin[S[idx].name]=margin
                killed=True; break
        if killed: continue
        axis=min(range(6),key=lambda i:depths[i])
        nd=list(depths); nd[axis]+=1; nd=tuple(nd)
        LS=[]; RS=[]
        for st in S:
            l,r=st.split(axis); LS.append(l); RS.append(r)
        stack.append((nd,RS)); stack.append((nd,LS))
        if visited>250000:
            raise RuntimeError('did not close')
    elapsed=time.time()-t0
    rho=(1-float(h0)**2)**0.5
    eps=float((Fraction(2,3)+delta)*(1-h0)/(1+h0))
    result={
        'status':'CERTIFIED_EMPTY',
        'delta':str(delta),'delta_float':float(delta),'X_delta':str(X),'X_delta_float':float(X),
        'h0':str(h0),'h0_float':float(h0),'rho_delta_lower_bound':rho,'epsilon_delta_lower_bound':eps,
        'outer_domain':{'p':['0','1/2'],'a':['-1','1'],'b':['-1','1'],'y':['0','1'],'u':['-1','1'],'v':['-1','1']},
        'statement':'No simultaneous solution of Ktilde=Cu=Cv=0, Cap<=0 and Nh>=0 exists on the outer box. Therefore no physical ordinary node with Xi_J<=X_delta and h>=h0 exists.',
        'visited_boxes':visited,'prune_counts':dict(prunes),'max_axis_depths':max_axis_depth,'deepest_depth_vector':list(deepest),
        'min_certified_Bernstein_sign_margin':{k:(None if q is None else str(q)) for k,q in minmargin.items()},
        'min_certified_Bernstein_sign_margin_float':{k:(None if q is None else float(q)) for k,q in minmargin.items()},
        'max_coefficient_error_bound_seen':{k:str(q) for k,q in maxerr_seen.items()},
        'max_coefficient_error_bound_seen_float':{k:float(q) for k,q in maxerr_seen.items()},
        'build':build,'subdivision_seconds':elapsed,
        'arithmetic_note':'Initial Bernstein coefficients are exact Fractions. Binary64 values are dyadic centers only. Their exact conversion errors and every propagated de Casteljau error radius are rational. At pruning, the selected floating minimum/maximum center is converted with Fraction.from_float and compared with the exact rational radius, so every sign decision is exact in Q. The child radius is E_parent + degree_axis * 2^-51 * Mbound.',
        'bernstein_note':'The polynomial range on a box lies in the convex hull of its exact Bernstein coefficients.',
        'hashes':{'Ktilde':expr_hash(Ktilde),'Cu':expr_hash(Cu),'Cv':expr_hash(Cv),'Nh':expr_hash(Nh),'Cap':expr_hash(Cap)}
    }
    with open(out_json,'w') as f: json.dump(result,f,indent=2,sort_keys=True)
    return result


def build_common():
    out=[]
    for name,E,sense in [('Ktilde',Ktilde,'zero'),('Cu',Cu,'zero'),('Cv',Cv,'zero')]:
        st,info=build_poly(name,E,sense); out.append((st,info)); print('built',name,info,flush=True)
    return out


if __name__=='__main__':
    common=build_common()
    for delta,h0,fn in [(Fraction(1,50),Fraction(24,25),'jnr_rational_interval_v3_delta002.json'),
                         (Fraction(1,20),Fraction(9,10),'jnr_rational_interval_v3_delta005.json')]:
        print('RUN',delta,h0,flush=True)
        R=run_certificate(delta,h0,fn,common)
        keys=['status','delta_float','h0_float','rho_delta_lower_bound','epsilon_delta_lower_bound','visited_boxes','prune_counts','max_axis_depths','min_certified_Bernstein_sign_margin_float','max_coefficient_error_bound_seen_float','subdivision_seconds']
        print(json.dumps({k:R[k] for k in keys},indent=2),flush=True)
