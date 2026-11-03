import numpy as np
from scipy.integrate import quad
from scipy.special import spherical_jn
def dT(s): # d/ds of T(s)=-2 j1(x)/x, x=s/sqrt3 ; = (2/sqrt3) j2(x)/x
    x=s/np.sqrt(3); return 2/np.sqrt(3)*spherical_jn(2,x)/x
def Kt_num(mu):
    re=quad(lambda s:2*dT(s)*np.cos(mu*s),1e-9,400,limit=2000)[0]
    im=quad(lambda s:2*dT(s)*np.sin(mu*s),1e-9,400,limit=2000)[0]
    return -1/3+re+1j*im
def Kt(mu):
    b=np.sqrt(3)*mu
    return 1-2*b*b-b*(1-b*b)*np.log(abs((1+b)/(1-b)))+1j*np.pi*b*(1-b*b)*(b<1)
for mu in (0.1,0.3,0.5,0.8):
    print(mu,Kt_num(mu),Kt(mu))
val=quad(lambda m:abs(Kt(m))**2,0,1,points=[1/np.sqrt(3)],limit=200)[0]
print('int |K|^2 =',val)
