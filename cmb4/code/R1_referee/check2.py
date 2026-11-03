import sympy as sp
x,eta,k,ep=sp.symbols('x eta k etap',positive=True)
X=k*eta/sp.sqrt(3)
TP=-2*(sp.sin(X)-X*sp.cos(X))/X**3
print('ODE',sp.simplify(sp.diff(TP,eta,2)+4/eta*sp.diff(TP,eta)+k**2/3*TP))
Td=-4*sp.cos(X)+8*sp.sin(X)/X+8*sp.cos(X)/X**2-8*sp.sin(X)/X**3
d00=-sp.Rational(2,3)*eta**2*k**2*TP-2*eta*sp.diff(TP,eta)-2*TP
print('delta from 00',sp.simplify(Td-d00))
TV=eta*(X**2*sp.sin(X)+2*X*sp.cos(X)-2*sp.sin(X))/X**3
V=-eta**2/2*(sp.diff(TP,eta)+TP/eta)
print('V',sp.simplify(TV-V))
# continuity: dT' = Phi' - (1/3) div v ; div v = -k^2 V
print('continuity',sp.simplify(sp.diff(Td/4,eta)-sp.diff(TP,eta)-sp.Rational(1,3)*k**2*TV))
print('Euler',sp.simplify(sp.diff(TV,eta)+TP+Td/4))
for f,nm in [(TP,'TP'),(Td,'Td'),(TV/eta,'TV/eta')]: print(nm,'->',sp.limit(f,eta,0))
# linear section forms
j0=lambda y: sp.sin(y)/y; j1=lambda y: sp.sin(y)/y**2-sp.cos(y)/y; j2=lambda y:(3/y**2-1)*sp.sin(y)/y-3*sp.cos(y)/y**2
print('TP=-2j1/x',sp.simplify(TP+2*j1(X)/X))
print('TdT',sp.simplify(Td/4-(-sp.cos(X)+2*j0(X)+TP)))
print('kTV',sp.simplify(k*TV-sp.sqrt(3)*(sp.sin(X)-2*j1(X))))
print('2dTP',sp.simplify(2*sp.diff(TP,eta)-4*k/sp.sqrt(3)*j2(X)/X))
# Green functions
xp=k*ep/sp.sqrt(3)
u1=lambda y:(sp.sin(y)-y*sp.cos(y))/y**3; u2=lambda y:(sp.cos(y)+y*sp.sin(y))/y**3
W=sp.simplify(u1(X)*sp.diff(u2(X),eta)-sp.diff(u1(X),eta)*u2(X))
print('W',W)
G=(u1(xp)*u2(X)-u1(X)*u2(xp))/W.subs(eta,ep)
G2=sp.sqrt(3)/k*xp/X**3*((1+X*xp)*sp.sin(X-xp)-(X-xp)*sp.cos(X-xp))
print('G forms',sp.simplify(sp.expand_trig(G-G2)))
L=lambda f: sp.diff(f,eta,2)+4/eta*sp.diff(f,eta)+k**2/3*f
print('LG',sp.simplify(L(G2)),'G(eta,eta)',sp.simplify(G2.subs(eta,ep)),'dG',sp.simplify(sp.diff(G2,eta).subs(eta,ep)))
print('G k->0',sp.simplify(sp.limit(G2,k,0)))
Gh=ep/eta*sp.sin(k*(eta-ep))/k
Lh=lambda f: sp.diff(f,eta,2)+2/eta*sp.diff(f,eta)+k**2*f
print('Gh',sp.simplify(Lh(Gh)),Gh.subs(eta,ep),sp.simplify(sp.diff(Gh,eta).subs(eta,ep)))
# normalisation algebra
e=sp.symbols('e')
m1,m2,m3,M1,M2,M3=sp.symbols('m1 m2 m3 M1 M2 M3') # M_i = sky monopole of m_i
M=1+e*m1+e**2/2*m2+e**3/6*m3; Mb=1+e*M1+e**2/2*M2+e**3/6*M3
th=sp.series(M/Mb,e,0,4).removeO()
for n in (1,2,3):
  c=sp.expand(sp.factorial(n)*th.coeff(e,n)); print('Theta',n,c)
L1,L2,L3=sp.symbols('L1 L2 L3')
ex=sp.series(sp.exp(e*L1+e**2/2*L2+e**3/6*L3),e,0,4).removeO()
print([sp.expand(sp.factorial(n)*ex.coeff(e,n)) for n in (1,2,3)])
# dilation members
a1,a2,c=sp.symbols('a1 a2 c')
A=sp.exp(e*a1+e**2/2*a2); B=sp.exp(e*c)
Psi=(A**4-1)/2; Phi=(1-A**2*B**2)/2; dl=A**-4-1
ser=lambda f,n: sp.expand(sp.factorial(n)*sp.series(f,e,0,3).removeO().coeff(e,n))
sol=sp.solve([ser(Psi,1)+sp.Rational(2,3)*c],[a1])
print('a1',sol)
s={a1:-c/3}
P2=ser(Phi,2).subs(s); S2=ser(Psi,2).subs(s); D2=ser(dl,2).subs(s)
a2s=sp.solve(sp.Eq(P2,S2),a2)[0]
print('Psi2=Phi2 member:',a2s,P2.subs(a2,a2s),D2.subs(a2,a2s), 'check delta2=-2Psi2+32/9c^2:',sp.simplify(D2+2*S2-sp.Rational(32,9)*c**2))
print('a2=0 member', P2.subs(a2,0),S2.subs(a2,0),D2.subs(a2,0))
