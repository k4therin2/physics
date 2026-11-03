import pickle, sympy as sp
ch=pickle.load(open('ch.pkl','rb'))
names={}
base={'TPhi':r'T_{\Phi}','TPhip':r'\dot T_{\Phi}','TdT':r'T_{\delta_T}','TdTp':r'\dot T_{\delta_T}','TV':r'T_{V}','TVp':r'\dot T_{V}',
 'Kf':r'\mathcal F','Dp':r'\mathcal D_{\perp}','Dz':r'\mathcal D_{r}','Xo':r'X_{o}','detao1':r'\delta\eta_{o}','detae1':r'\delta\eta_{e}'}
for b,l in base.items():
    for a in (1,2):
        if l.endswith('}'): names[sp.Symbol('%s_%d'%(b,a))]=l[:-1]+',%d}'%a
        else: names[sp.Symbol('%s_%d'%(b,a))]=l+'_{%d}'%a
names[sp.Symbol('E1')]='E_1';names[sp.Symbol('E2')]='E_2';names[sp.Symbol('EK')]='E_K'
names[sp.Symbol('Pperp')]=r'P_\perp'
for n,l in {'Phi2':r'\Phi_2','Phi2p':r'\dot\Phi_2','Psi2':r'\Psi_2','Psi2p':r'\dot\Psi_2','dT2':r'\delta_{T2}','V2':'V_2',
 'KvT12':r'K_{v^T}^{(12)}','KvT21':r'K_{v^T}^{(21)}','Kom12':r'K_\omega^{(12)}','Kom21':r'K_\omega^{(21)}','H':r'\mathcal H','Hp':r'\dot{\mathcal H}'}.items():
    names[sp.Symbol(n)]=l
for n,l in {'k1':'k_1','k2':'k_2','mu':r'\mu','mu1':r'\mu_1','mu2':r'\mu_2','K':'K','Kn':'K_n'}.items():
    names[sp.Symbol(n,real=True)]=l
title={'[ray] local | fields:2':('ray','second-order fields (C: $\\Phi_2,\\Psi_2$ ISW, $\\omega_2$, $h_2$)'),
'[ray] lensing (transverse Dx1) | fields:1':('ray','lensing, transverse $\\Delta\\vx_1$ (D)'),
'[ray] time delay (radial Dx1) | fields:1':('ray','time delay, radial $\\Delta\\vx_1$ (D)'),
'[ray] transverse direction k1 (aberration+deflection) | fields:1':('ray','aberration + deflection, transverse $\\vk_1$ (D)'),
'[ray] local | fields:1+1':('ray','local $1\\times1$ (C)'),
'[obs] observer clock | fields:1':('obs','observer clock (D)'),
'[emi] emission-surface shift | fields:1':('emi','emission-surface shift (D)'),
'[obs] local | fields:1+1':('obs','local $1\\times1$ (B, endpoint)'),
'[obs] local | fields:2':('obs','second-order fields (C: $\\Psi_2,\\delta_{T2},V_2,v^T_2+\\omega_2$)'),
'[obs] observer displacement | fields:1':('obs','observer displacement (D)'),
'[emi] local | fields:1+1':('emi','local $1\\times1$ (B, endpoint)'),
'[emi] local | fields:2':('emi','second-order fields (C: $\\Psi_2,\\delta_{T2},V_2,v^T_2+\\omega_2$)'),
'[emi] lensing (transverse Dx1) | fields:1':('emi','lensing at emission (D)'),
'[emi] time delay (radial Dx1) | fields:1':('emi','time delay at emission (D)'),
'[emi] transverse direction k1 (aberration+deflection) | fields:1':('emi','aberration + deflection at emission (D)')}
where={'ray':r'$\int_{\eta_*}^{\eta_o}\!\dd\eta$ of','obs':r'at $\eta=\eta_o$, $\chi=0$ ($E_a=1$):','emi':r'at $\eta=\eta_*$, $\chi=\chi_*$:'}
out=[]
for k,(w,e) in ch.items():
    t=title[k]
    out.append(r'\paragraph{[%s] %s.} %s'%(t[0],t[1],where[t[0]]))
    out.append(r'\begin{dmath*}'+'\n'+sp.latex(sp.expand(e),symbol_names=names,order='none')+'\n'+r'\end{dmath*}')
open('k2body.tex','w').write('\n'.join(out)+'\n')
print(open('k2body.tex').read()[:1500])
