import sympy as sp
e=sp.symbols('epsilon')
# 1. zeta expansions
d1,d2,d3,P1,P2,P3=sp.symbols('d1 d2 d3 P1 P2 P3')
d=e*d1+e**2/2*d2+e**3/6*d3; P=e*P1+e**2/2*P2+e**3/6*P3
z=sp.Rational(1,4)*sp.log(1+d)+sp.Rational(1,2)*sp.log(1-2*P)
s=sp.series(z,e,0,4).removeO()
for n in (1,2,3): print('zeta',n, sp.expand(sp.factorial(n)*s.coeff(e,n)))
# 2. K_Q forms, K_Phi2 forms
k,q,p,mu=sp.symbols('k q p mu',positive=True)
psub={p:sp.sqrt(k**2-2*k*q*mu+q**2)}
KQa=-(13*k**4-6*k**2*(p**2+q**2)+9*(p**2-q**2)**2)/(9*k**4)
KQb=-sp.Rational(4,9)*(4*k**2-6*k*q*mu+9*q**2*mu**2-3*q**2)/k**2
print('KQ forms',sp.simplify(KQa.subs(psub)-KQb))
KPa=-(7*k**4-2*k**2*(p**2+q**2)+3*(p**2-q**2)**2)/(9*k**4)
KPb=-sp.Rational(4,9)*(2*k**2-2*k*q*mu+3*q**2*mu**2-q**2)/k**2
print('KPhi forms',sp.simplify(KPa.subs(psub)-KPb), sp.simplify(KQa/3-sp.Rational(8,27)-KPa))
# position space Phi2(0): Phi1=-2/3 zeta -> kernel per zeta zeta; qv, pv vectors
qx,qz,px,pz=sp.symbols('qx qz px pz',real=True)
qv=sp.Matrix([qx,0,qz]); pv=sp.Matrix([px,0,pz]); kv=qv+pv
dot=lambda a,b:(a.T*b)[0]
c=sp.Rational(4,9) # Phi1 Phi1 kernel
ksq=dot(kv,kv)
def sym(f): return (f+f.subs({qx:px,px:qx,qz:pz,pz:qz},simultaneous=True))/2
Phi2pos = -2*c + (-1)*(-1/ksq)*(-dot(qv,pv))*c + 3*(1/ksq**2)*sym(dot(kv,qv)*dot(kv,pv))*c
# careful: d_i d_j (d_i Phi d_j Phi) -> (i k_i)(i k_j)(i q_i)(i p_j) = k.q k.p
Phi2pos = -2*c + (-1)*(-1/ksq)*(-dot(qv,pv))*c + 3*(1/ksq**2)*(dot(kv,qv)*dot(kv,pv))*c
inv={k:sp.sqrt(ksq),q:sp.sqrt(dot(qv,qv)),p:sp.sqrt(dot(pv,pv))}
print('Phi2 pos vs kernel',sp.simplify(Phi2pos-KPa.subs(inv)))
Psi2pos= 4*c + 3*(-1/ksq)*(-dot(qv,pv))*c - 9*(1/ksq**2)*(dot(kv,qv)*dot(kv,pv))*c
print('Psi2-Phi2 pos vs -KQ',sp.simplify(Psi2pos+KQa.subs(inv)))
# 3. C3 forms
C3a=-16*(2*k**5+12*k**4*q*mu-42*k**3*q**2*mu**2-10*k**3*q**2+36*k**2*q**3*mu**3+32*k**2*q**3*mu+3*k*q**4*mu**4-42*k*q**4*mu**2+k*q**4+8*q**5*mu)/(9*k*(k**2-2*k*q*mu+q**2)**2)
C3b=-(35*k**8+12*k**6*p**2-84*k**6*q**2+66*k**4*p**4-52*k**4*p**2*q**2+66*k**4*q**4-84*k**2*p**6+84*k**2*p**4*q**2+20*k**2*p**2*q**4-20*k**2*q**6+3*p**8-12*p**6*q**2+18*p**4*q**4-12*p**2*q**6+3*q**8)/(9*k**4*p**4)
print('C3 forms',sp.simplify(C3a-C3b.subs(psub)))
print('C3 q->0',sp.limit(C3a,q,0))
num=sp.numer(sp.together(C3a))
print('C3 numerator at q=k,mu=1:',sp.factor(num.subs({q:k,mu:1})))
t=sp.symbols('t')
print('C3 near q=k, mu=1-t:',sp.series(C3a.subs({q:k,mu:1-t}),t,0,1))
print('C3 at q=k exact in mu:',sp.factor(sp.simplify(C3a.subs(q:=k,k) if False else C3a.subs({q:k}))))
# muaverage of C3 at q=k? integrable?
