import sympy as sp
e=sp.symbols('e')
P1,P2,F1,F2=sp.symbols('Psi1 Psi2 Phi1 Phi2')
v1=sp.symbols('v1x v1y v1z'); v2=sp.symbols('v2x v2y v2z'); w2=sp.symbols('w2x w2y w2z')
hxx,hyy,hxy,hxz,hyz=sp.symbols('hxx hyy hxy hxz hyz'); hzz=-hxx-hyy
H=sp.Matrix([[hxx,hxy,hxz],[hxy,hyy,hyz],[hxz,hyz,hzz]])
Psi=e*P1+e**2/2*P2; Phi=e*F1+e**2/2*F2
v=[e*v1[i]+e**2/2*v2[i] for i in range(3)]; w=[e**2/2*w2[i] for i in range(3)]
g=sp.zeros(4,4); g[0,0]=-(1+2*Psi)
for i in range(3):
  g[0,i+1]=g[i+1,0]=w[i]
  for j in range(3): g[i+1,j+1]=(1-2*Phi)*(1 if i==j else 0)+e**2/2*H[i,j]
O=3
tr=lambda f: sp.series(sp.expand(f),e,0,O).removeO()
U=sp.Matrix([1]+v)
norm=tr((U.T*g*U)[0])
u0=tr(1/sp.sqrt(-norm))
u=(u0*U).applyfunc(tr)
dot=lambda a,b: tr((a.T*g*b)[0])
basis=[]
for i in range(3):
  d=sp.zeros(4,1); d[i+1]=1
  vec=(d+dot(u,d)*u).applyfunc(tr)
  for b in basis: vec=(vec-dot(b,vec)*b).applyfunc(tr)
  nrm=tr(1/sp.sqrt(dot(vec,vec)))
  basis.append((vec*nrm).applyfunc(tr))
E3=basis[2]
# triad check vs doc for a=3 (z) and a=1
print('e_3^0 doc diff:',sp.simplify(E3[0]-(e*v1[2]+e**2*(sp.Rational(1,2)*(v2[2]+w2[2])-(F1+2*P1)*v1[2]))))
docspatial=1+e*F1+e**2*(sp.Rational(3,2)*F1**2+F2/2-hzz/4+v1[2]**2/2)
print('e_3^z diff:',sp.simplify(E3[3]-docspatial))
print('e_3^x (i<a):',sp.simplify(E3[1]), ' doc: -1/2 h_{2,zx} e^2 (with h=e^2/2 h2 => coefficient check)')
E1=basis[0]
print('e_1^z (i>a):',sp.simplify(E1[3]),' doc: e^2 v1x v1z')
kv=[tr(sp.series((u[i+1]-E3[i+1])/(u[0]-E3[0]),e,0,O).removeO()) for i in range(3)]
print('k^x 1st:',sp.expand(kv[0]).coeff(e,1),' 2nd/ (eps^2 coeff):',sp.expand(kv[0]).coeff(e,2))
print('k^z 1st:',sp.expand(kv[2]).coeff(e,1))
kz2=sp.expand(kv[2]).coeff(e,2)
doc=-P2/2-F2/2-sp.Rational(3,2)*F1**2-F1*P1+P1**2/2+hzz/4-w2[2]/2+(v1[0]**2+v1[1]**2)/2
print('k^z eps^2 minus doc:',sp.simplify(kz2-doc))
kx2=sp.expand(kv[0]).coeff(e,2)
docx=v2[0]/2+v1[0]*v1[2]+hxz/2
print('k^x eps^2 minus doc:',sp.simplify(kx2-docx))
