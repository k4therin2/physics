import sympy as sp
k,q,mu,eta=sp.symbols('k q mu eta',positive=True)
Tq,Tp,dTq,dTp=sp.symbols('Tq Tp dTq dTp')
p2=k**2-2*k*q*mu+q**2
app=-(Tp*Tq*(13*k**4-6*k**2*p2-6*k**2*q**2+9*p2**2-18*p2*q**2+9*q**4))/(4*k**4) \
 -(eta*(Tp*dTq+Tq*dTp)+eta**2*dTp*dTq)*(-k**4-2*k**2*p2-2*k**2*q**2+3*p2**2-6*p2*q**2+3*q**4)/(4*k**4)
main=Tq*Tp*(-4+6*mu*q/k+3*q**2/k**2-9*mu**2*q**2/k**2)+(eta*(Tp*dTq+Tq*dTp)+eta**2*dTq*dTp)*(2*mu*q/k+q**2/k**2-3*mu**2*q**2/k**2)
print(sp.simplify(sp.expand(app-main)))
# K_Q(0) from app with T->-2/3, eta dT ->0
print(sp.simplify(app.subs({Tq:-sp.Rational(2,3),Tp:-sp.Rational(2,3),dTq:0,dTp:0})+sp.Rational(4,9)*(4*k**2-6*k*q*mu+9*q**2*mu**2-3*q**2)/k**2))
