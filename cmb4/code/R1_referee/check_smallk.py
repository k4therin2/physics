import sympy as sp
k,eta,es,eo=sp.symbols('k eta eta_s eta_o',positive=True)
N=8
X=lambda h:k*h/sp.sqrt(3)
def ser(f): return sp.series(f,k,0,N).removeO()
j0=lambda y: sp.sin(y)/y; j1=lambda y: sp.sin(y)/y**2-sp.cos(y)/y
TP=lambda h:-2*j1(X(h))/X(h)
A=lambda h: -sp.cos(X(h))+2*j0(X(h))-4*j1(X(h))/X(h)
kTV=lambda h: sp.sqrt(3)*(sp.sin(X(h))-2*j1(X(h)))
w=k*(eo-es)
integrand=ser(2*sp.diff(TP(eta),eta)*j0(k*(eo-eta)))
I0=sp.integrate(sp.expand(integrand),(eta,es,eo))
D0=sp.expand(ser(A(es)*j0(w)+kTV(es)*j1(w)-A(eo))+I0)
print('Delta0 series:',sp.factor(sp.collect(D0,k).coeff(k,2)), '| k^4:',sp.factor(D0.coeff(k,4)))
# Delta1: A j1(w) - kTV j1'(w) + 1/3 kTV(eo) + I1
j1p=lambda y: sp.diff(j1(sp.Symbol('u')),sp.Symbol('u')).subs(sp.Symbol('u'),y)
I1=sp.integrate(sp.expand(ser(2*sp.diff(TP(eta),eta)*j1(k*(eo-eta)))),(eta,es,eo))
D1=sp.expand(ser(A(es)*j1(w)-kTV(es)*j1p(w)+kTV(eo)/3)+I1)
print('Delta1 k^1:',sp.factor(D1.coeff(k,1)),'| k^3:',sp.factor(D1.coeff(k,3)))
