import numpy as np
from scipy.special import sph_harm_y, spherical_jn
sph_harm=lambda m,l,ph,th: sph_harm_y(l,m,th,ph)
from scipy.integrate import dblquad
import math
x=2.7
def coef(l,lam):
    def f(th,ph,part):
        e=np.array([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)])
        e1=np.array([1,1j,0])/np.sqrt(2)
        Q=(e@e1) if lam==1 else (e@e1)**2
        val=np.conj(sph_harm(lam,l,ph,th))*Q*np.exp(-1j*x*np.cos(th))*np.sin(th)
        return val.real if part==0 else val.imag
    re=dblquad(lambda th,ph:f(th,ph,0),0,2*np.pi,0,np.pi)[0]
    im=dblquad(lambda th,ph:f(th,ph,1),0,2*np.pi,0,np.pi)[0]
    return re+1j*im
for lam in (1,2):
  for l in (2,3,4):
    c=coef(l,lam)
    j=spherical_jn(l,x)
    if lam==1: D=-1j*np.sqrt(l*(l+1)/2)*j/x
    else: D=-0.5*np.sqrt(math.factorial(l+2)/math.factorial(l-2))*j/x**2
    pred=(-1j)**l*np.sqrt(4*np.pi*(2*l+1))*D
    print(lam,l,c,pred)
