# Independent check of the TT source normalisation K_h (sec_einstein eq second_kernels)
import sympy as sp
t,x,y,z=sp.symbols('eta x y z',positive=True)
e=sp.symbols('epsilon')
qx,qz,pz=sp.symbols('qx qz pz',positive=True)
A=sp.Function('A')(t); B=sp.Function('B')(t)
Ea=sp.exp(sp.I*(qx*x+qz*z)); Eb=sp.exp(sp.I*(-qx*x+pz*z))
q2=qx**2+qz**2; p2=qx**2+pz**2
# first-order relations (per unit zeta) : Phi=Psi=A Ea + B Eb ; V=-eta^2/2(Phi'+Phi/eta); delta=-2/3 eta^2 k^2 Phi -2 eta Phi' -2 Phi
def V(F): return -t**2/2*(sp.diff(F,t)+F/t)
def dl(F,kk2): return -sp.Rational(2,3)*t**2*kk2*F-2*t*sp.diff(F,t)-2*F
Phi=e*(A*Ea+B*Eb)
Vf=e*(V(A)*Ea+V(B)*Eb)
delta=e*(dl(A,q2)*Ea+dl(B,p2)*Eb)
X=[t,x,y,z]
v=[sp.diff(Vf,xx) for xx in (x,y,z)]
gh=sp.diag(-(1+2*Phi),1-2*Phi,1-2*Phi,1-2*Phi)
g=t**2*gh
def trunc(f): return sp.series(sp.expand(f),e,0,3).removeO()
ginv=sp.diag(*[trunc(1/g[i,i]) for i in range(4)])
# u^mu
u0h=sp.symbols('u0h')
vv=sum(vi**2 for vi in v)
u0hat=trunc(1/sp.sqrt(1+2*Phi-(1-2*Phi)*vv))
uup=[u0hat/t]+[u0hat*vi/t for vi in v]
udn=[trunc(sum(g[m,n]*uup[n] for n in range(4))) for m in range(4)]
# Christoffel
Gam=[[[0]*4 for _ in range(4)] for _ in range(4)]
for a in range(4):
  for b in range(4):
    for c in range(b,4):
      val=sum(ginv[a,d]*(sp.diff(g[d,b],X[c])+sp.diff(g[d,c],X[b])-sp.diff(g[b,c],X[d])) for d in range(4))/2
      Gam[a][b][c]=Gam[a][c][b]=trunc(val)
def Ric(b,c):
  r=0
  for a in range(4):
    r+=sp.diff(Gam[a][b][c],X[a])-sp.diff(Gam[a][b][a],X[c])
    for d in range(4):
      r+=Gam[a][a][d]*Gam[d][b][c]-Gam[a][c][d]*Gam[d][b][a]
  return trunc(r)
R={}
for (b,c) in [(1,1),(2,2)]: R[(b,c)]=Ric(b,c)
# Ricci scalar needed: G_ij = R_ij - 1/2 g_ij R ; for xx-yy combination g_xx=g_yy so R scalar drops out!
Exx_m_yy = R[(1,1)]-R[(2,2)] - (3/t**4)*(1+delta)*sp.Rational(4,3)*(udn[1]**2-udn[2]**2)  # g_xx-g_yy=0
c2=sp.expand(Exx_m_yy).coeff(e,2)
c2=sp.expand(c2)
cross=sp.simplify(sp.expand(c2*Ea**-1*Eb**-1).subs({x:0,y:0,z:0}))  # crude: need only the Ea*Eb terms
# extract properly: collect terms proportional to Ea*Eb
terms=sp.Add.make_args(c2)
cc=0
for tm in terms:
  r=sp.simplify(tm/(Ea*Eb))
  if not r.has(x) and not r.has(z): cc+=r
cc=sp.simplify(cc)
# impose ODEs for A,B: A''=-4/eta A' - q^2/3 A
sub={sp.Derivative(A,(t,2)):-4/t*sp.diff(A,t)-q2/3*A, sp.Derivative(B,(t,2)):-4/t*sp.diff(B,t)-p2/3*B}
def red(f):
  for _ in range(3):
    f=f.subs({sp.Derivative(A,(t,3)):sp.diff(-4/t*sp.diff(A,t)-q2/3*A,t),sp.Derivative(B,(t,3)):sp.diff(-4/t*sp.diff(B,t)-p2/3*B,t)}).subs(sub)
  return sp.simplify(sp.expand(f))
cc=red(cc)
Kh_doc=4*(2*A*B+(t*sp.diff(A,t)+A)*(t*sp.diff(B,t)+B))
pred=-2*cc/qx**2
print('ratio  (-2 c / qx^2) / K_h_doc =', sp.simplify(pred/Kh_doc))
print('difference', sp.simplify(sp.expand(pred-Kh_doc)))
