import sympy as sp
k,q,mu,p,eta,x,eo,K,c,ks=sp.symbols('k q mu p eta x eta_o K c k_s',positive=True)
# KQ0 and KPhi20 forms
p2=k**2-2*k*q*mu+q**2
KQ_a=-(13*k**4-6*k**2*(p2+q**2)+9*(p2-q**2)**2)/(9*k**4)
KQ_b=-sp.Rational(4,9)*(4*k**2-6*k*q*mu+9*q**2*mu**2-3*q**2)/k**2
print('KQ forms',sp.simplify(KQ_a-KQ_b))
KP_a=-(7*k**4-2*k**2*(p2+q**2)+3*(p2-q**2)**2)/(9*k**4)
KP_b=-sp.Rational(4,9)*(2*k**2-2*k*q*mu+3*q**2*mu**2-q**2)/k**2
print('KPhi forms',sp.simplify(KP_a-KP_b), sp.simplify(KQ_b/3-sp.Rational(8,27)-KP_b))
print('KPhi smallk', sp.expand(KP_b))
# zeta2(0)=0 check: delta2 = 32/9 - 2 Psi2, Psi2=Phi2-KQ
Phi2=sp.Symbol('Phi2')
d2=sp.Rational(32,9)-2*(Phi2-KQ_b); d1=sp.Rational(4,3); P1=-sp.Rational(2,3)
z2=d2/4-d1**2/4-Phi2-2*P1**2
print('Phi2 solve',sp.simplify(sp.solve(z2,Phi2)[0]-KP_b))
# C3
C3=-16*(2*k**5+12*k**4*q*mu-42*k**3*q**2*mu**2-10*k**3*q**2+36*k**2*q**3*mu**3+32*k**2*q**3*mu+3*k*q**4*mu**4-42*k*q**4*mu**2+k*q**4+8*q**5*mu)/(9*k*p2**2)
C3k=sp.factor(C3.subs(q,k)); print('C3(q=k)',C3k, sp.simplify(C3k+4*(3*mu**3+39*mu**2-45*mu+7)/(9*(mu-1))))
e=sp.symbols('e',positive=True)
print('C3 pole series', sp.series(C3k.subs(mu,1-e),e,0,1))
print('C3 q->0', sp.limit(C3,q,0))
print('6*1/2*(4/3)^2*(-2/3)=',6*sp.Rational(1,2)*sp.Rational(16,9)*(-sp.Rational(2,3)))
# clock and displacement
xo=sp.symbols('x_o',positive=True)
TP=-2*(sp.sin(x)-x*sp.cos(x))/x**3
TV_over_eta=(x**2*sp.sin(x)+2*x*sp.cos(x)-2*sp.sin(x))/x**3
# eta = sqrt3 x/k
I1=sp.integrate(x*TP,(x,0,xo))*3/k**2
deta=sp.simplify(-I1/(sp.sqrt(3)*xo/k))
print('deta_o', sp.simplify(deta-6/(k**2*(sp.sqrt(3)*xo/k))*(1-sp.sin(xo)/xo)))
Xo=sp.integrate(x*TV_over_eta,(x,0,xo))*3/k**2
print('X_o', sp.simplify(Xo-3/k**2*(2*sp.sin(xo)/xo-sp.cos(xo)-1)))
# Theta3 normalisation
eps=sp.symbols('epsilon'); m1,m2,m3,a1,a2,a3=sp.symbols('m1 m2 m3 a1 a2 a3')
M=1+eps*m1+eps**2/2*m2+eps**3/6*m3; Mb=1+eps*a1+eps**2/2*a2+eps**3/6*a3
s=sp.series(M/Mb,eps,0,4).removeO()
th3=sp.expand(6*s.coeff(eps,3)); th2=sp.expand(2*s.coeff(eps,2))
print('Th2',th2); print('Th3',th3)
# m_n from exp
L1,L2,L3=sp.symbols('L1 L2 L3')
ss=sp.series(sp.exp(eps*L1+eps**2/2*L2+eps**3/6*L3),eps,0,4).removeO()
print('m2',sp.expand(2*ss.coeff(eps,2)),' m3',sp.expand(6*ss.coeff(eps,3)))
