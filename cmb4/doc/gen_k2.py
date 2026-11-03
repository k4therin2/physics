import pickle, sympy as sp, sys
from sympy import *
P='/home/user/physics/cmb4/code/A2_second_order/symbolic_kernels.pkl'
d=pickle.load(open(P,'rb'))
ch={k:(w,sp.sympify(eval(e,sp.__dict__))) for k,(w,e) in d['channels'].items()}
S=sp.Symbol
k1,k2,mu,mu1,mu2,K,Kn=[S(n,real=True) for n in 'k1 k2 mu mu1 mu2 K Kn'.split()]
EK=S('EK');Pp=S('Pperp')
sec={n:S(n) for n in 'Phi2 Phi2p Psi2 Psi2p dT2 V2 KvT12 KvT21 Kom12 Kom21 H Hp'.split()}
def nPk(s):
    km={1:k1,2:k2}[s]; mm={1:mu1,2:mu2}[s]
    return km*mm-Kn*(km**2+k1*k2*mu)/K**2
nom=(nPk(1)*sec['Kom12']+nPk(2)*sec['Kom21'])/2
nvT=(nPk(1)*sec['KvT12']+nPk(2)*sec['KvT21'])/2
hnn=lambda HH:(nPk(1)**2-Rational(1,2)*(1-Kn**2/K**2)*(k1**2-(k1**2+k1*k2*mu)**2/K**2))*HH
# claimed decompositions
ray2=EK*(sec['Phi2p']-sec['Psi2p']+2*I*Kn*sec['Psi2']-hnn(sec['Hp'])/2+I*Kn*nom)
obs2=EK*(sec['Psi2']-sec['dT2']+I*Kn*sec['V2']+nvT+nom)
def chk(name,claim):
    w,e=ch[name]; r=sp.simplify(sp.expand(e-claim)); print(name,'residual:',r); assert r==0
# map symbol names
fix=lambda e:e
chk('[ray] local | fields:2',ray2)
chk('[obs] local | fields:2',obs2)
chk('[emi] local | fields:2',-obs2)
# also check first-order channels as quoted in text
def s(n,a): return S('%s_%d'%(n,a))
E={1:S('E1'),2:S('E2')}
def sym12(f): return f(1,2)+f(2,1)
km={1:k1,2:k2}; mm={1:mu1,2:mu2}
chk('[ray] lensing (transverse Dx1) | fields:1', sym12(lambda a,b:-2*I*Pp*s('Dp',a)*s('TPhi',b)*km[b]*mm[b]*E[b]))
chk('[ray] time delay (radial Dx1) | fields:1', sym12(lambda a,b:-2*s('Dz',a)*s('TPhi',b)*km[b]**2*mm[b]**2*E[b]))
chk('[ray] transverse direction k1 (aberration+deflection) | fields:1', sym12(lambda a,b:2*Pp*s('Kf',a)*s('TPhi',b)*E[b]))
chk('[ray] local | fields:1+1', 4*E[1]*E[2]*sym12(lambda a,b:s('TPhi',a)*s('TPhip',b)))
U=lambda a,b: s('TPhi',a)+I*km[a]*mm[a]*s('TV',a)-s('TdT',a)   # Psi + n.v - dT  (n.v1 = i k mu T_V)
chk('[emi] lensing (transverse Dx1) | fields:1', sym12(lambda a,b:Pp*s('Dp',a)*E[b]*U(b,a)))
chk('[emi] time delay (radial Dx1) | fields:1', sym12(lambda a,b:-I*km[b]*mm[b]*s('Dz',a)*E[b]*U(b,a)))
chk('[emi] transverse direction k1 (aberration+deflection) | fields:1', sym12(lambda a,b:-Pp*s('Kf',a)*s('TV',b)*E[b]))
# emission shift: claim -detae1_a * E_b * [ D U_b ] with D = d_eta - n.grad -> dot - i k mu
DU=lambda b: (s('TPhip',b)+I*km[b]*mm[b]*s('TVp',b)-s('TdTp',b)) - I*km[b]*mm[b]*U(b,None)
J=lambda b:2*I*km[b]*mm[b]*s('TPhi',b)
chk('[emi] emission-surface shift | fields:1', sym12(lambda a,b:-s('detae1',a)*E[b]*(DU(b)+J(b))))
chk('[obs] observer clock | fields:1', sym12(lambda a,b:s('detao1',a)*E[b]*(DU(b)+J(b)+I*km[b]*mm[b]*U(b,None))))
chk('[obs] observer displacement | fields:1', sym12(lambda a,b:-s('Xo',a)*k1*k2*mu*E[b]*U(b,None)).subs(k1*k2*mu,Pp+k1*k2*mu1*mu2))
w,e=ch['[obs] local | fields:1+1']; print('obs11', e)
w,e=ch['[emi] local | fields:1+1']; print('emi11', e)
pickle.dump(ch,open('ch.pkl','wb'))
