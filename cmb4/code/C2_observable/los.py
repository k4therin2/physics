"""Milestone 3: S = int_{eta_e}^{eta_o} deta sigma_kk/(-u.k) on the perturbed ray between perturbed endpoints,
expanded about the background ray, through eps^3.  Lambda = -S = eps L1 + eps^2 L2/2 + eps^3 L3/6.

Representation. Every integrand/endpoint expression is a sympy expression in
  * fields Psi_n, Phi_n, dT_n, v_i_n, om_i_n, h_ij_n as Functions of (eta,x,y,z)  -> to be evaluated at
      where='ray': (eta, xbar(eta) = xbar_o + n(etabar_o - eta)), integrated  int_{eta_*}^{etabar_o} deta
      where='obs': (etabar_o, xbar_o)          where='emi': (eta_*, xbar_*)
  * ray functionals k1_i(eta), k2_i(eta), Dx1_i(eta), Dx2_i(eta)  (definitions: ray.definitions; at 'obs' they are
      replaced by their observer values; at 'emi' they appear as k1_i(eta_star) etc.)
  * observer-event shifts deta_o1, deta_o2, dxo1_i, dxo2_i (symbols, direction independent)
  * emission-time shifts deta_e1, deta_e2 (symbols) with definitions in `emission` (functions of n through xbar_*)
  * n_i (sky direction; |n|=1 is NOT imposed in the general-n expressions).
The 'split' display version sets n = z-hat so that k1_3 = -(Phi1+Psi1) (null condition), Dx1_3 = radial (time delay),
k1_{1,2}, Dx1_{1,2} = transverse (aberration/deflection, lensing).
"""
import sympy as sp, pickle, time, itertools
from ser import Ser, N, eps
from fields import build, X, eta, x, y, z, F, fname, dT, v
import ray
from ray import n, kfun, Dxfun, G1, G2half, k1_init, k2_init_half, deta_o1, deta_o2, dxo1, dxo2, unit_rules

t0 = time.time()
d = build(verbose=False)
g, Gam, fl, rules = d['g'], d['Gam'], d['fl'], d['rules']
sig, ul, theta = fl['sig'], fl['ul'], fl['theta']
kser = ray.kser
eta_s = sp.Symbol('eta_star'); eta_o = sp.Symbol('eta_o_bar')
deta_e1, deta_e2 = sp.symbols('deta_e1 deta_e2')

def contract2(Tm):
    return sum((Tm[m][q]*kser[m]*kser[q] for m in range(4) for q in range(4) if not Tm[m][q].iszero()), Ser.const(0))
muk = -sum((ul[m]*kser[m] for m in range(4)), Ser.const(0))     # -u.k
I_sigma = contract2(sig)/muk                                      # sigma_kk/(-u.k)   [shear form]
Gam0kk = contract2(Gam[0])
lnT = (1 + dT).log()
B_direct = muk.log() - lnT                                        # boundary function of route (i): ln(-u.k) - ln That
I_direct = -Gam0kk                                                # integrand of route (i)
print('integrands built, t=%.1f' % (time.time()-t0))

# ---------- positional Taylor expansion about the background ray ----------
Dx = [Ser([0, Dxfun[1][i], Dxfun[2][i]/2, 0]) for i in range(3)]
def pos_expand(Q):
    """Q(x + Dx) = sum_m (Dx.d)^m Q / m!  (Dx = O(eps))"""
    out = Q; term = Q
    for m in range(1, N+1):
        term = sum((Dx[i]*term.diff(X[i+1]) for i in range(3)), Ser.const(0))*sp.Rational(1, m)
        if term.iszero(): break
        out = out + term
    return out

# ---------- total derivative along the ray in the background-ray representation ----------
G2tot_half = [G2half[i] + sum(Dxfun[1][j]*sp.diff(G1[i], X[j+1]) for j in range(3)) for i in range(3)]
Drep = {}
for i in range(3):
    Drep[sp.Derivative(kfun[1][i], eta)] = G1[i]
    Drep[sp.Derivative(kfun[2][i], eta)] = 2*G2tot_half[i]
    Drep[sp.Derivative(Dxfun[1][i], eta)] = kfun[1][i]
    Drep[sp.Derivative(Dxfun[2][i], eta)] = kfun[2][i]
def Dtot_expr(e):
    r = sp.diff(e, eta) - sum(n[i]*sp.diff(e, X[i+1]) for i in range(3))
    return sp.expand(r.xreplace(Drep))
def Dtot(Q, upto=N): return Ser([Dtot_expr(t) if i <= upto else sp.Integer(0) for i, t in enumerate(Q.c)])
def trunc(Q, upto): return Ser([t if i <= upto else sp.Integer(0) for i, t in enumerate(Q.c)])

def shift_eta(Q, dser):
    """R(eta + deta) = sum_m deta^m/m! Dtot^m R, deta = dser (Ser, O(eps))."""
    out = Q; term = Q; p = Ser.const(1)
    for m in range(1, N+1):
        term = Dtot(term, upto=N-m); p = p*dser
        add = term*p*sp.Rational(1, sp.factorial(m))
        if add.iszero(): break
        out = out + add
    return out

# ---------- endpoint substitutions ----------
sub_obs_fun = {}   # functionals at etabar_o
for i in range(3):
    sub_obs_fun[kfun[1][i]] = k1_init[i]
    sub_obs_fun[Dxfun[1][i]] = dxo1[i] + n[i]*deta_o1
k2o_half = [k2_init_half[i] + deta_o1*sp.diff(k1_init[i], eta) + sum(dxo1[j]*sp.diff(k1_init[i], X[j+1]) for j in range(3)) - deta_o1*G1[i] for i in range(3)]
for i in range(3):
    sub_obs_fun[kfun[2][i]] = 2*k2o_half[i]
    sub_obs_fun[Dxfun[2][i]] = 2*(dxo2[i]/2 + n[i]*deta_o2/2 - deta_o1*k1_init[i])
sub_emi_fun = {}
for i in range(3):
    for f_ in (kfun[1][i], kfun[2][i], Dxfun[1][i], Dxfun[2][i]):
        sub_emi_fun[f_] = f_.subs(eta, eta_s)
def at_obs(e): return sp.expand(sp.expand(e).xreplace(sub_obs_fun))
def at_emi(e): return sp.expand(sp.expand(e).xreplace(sub_emi_fun))

# ---------- emission time shift: 1 + dT(eta_e, x(eta_e)) = eta_e/eta_* ----------
RT = pos_expand(dT)
dser_e = Ser([0, deta_e1, deta_e2/2, 0])
lhs = shift_eta(RT, dser_e)            # dT on the perturbed ray at eta_* + deta_e  (fields at 'emi')
rhs = dser_e*(1/eta_s)
emission = {}
e1 = sp.expand(at_emi(lhs.c[1] - rhs.c[1]))
emission[deta_e1] = sp.solve(e1, deta_e1)[0]
e2 = sp.expand(at_emi(lhs.c[2] - rhs.c[2]).subs(deta_e1, emission[deta_e1]))
emission[deta_e2] = sp.expand(sp.solve(e2, deta_e2)[0])       # deta_e2 itself (dser_e carries deta_e2/2)
print('emission shifts: deta_e1 =', emission[deta_e1]); print('   deta_e2 =', emission[deta_e2]); print('t=%.1f' % (time.time()-t0))

# ---------- assemble S for a given integrand (and optional boundary function) ----------
dser_o = Ser([0, deta_o1, deta_o2/2, 0])
def assemble(I, B=None):
    """S = int_{eta_e}^{eta_o} I deta + [B]_e^o on the perturbed ray. Returns dict order-> list of (where, expr)."""
    Fq = pos_expand(I); print('   pos_expand done, sizes', [len(sp.Add.make_args(c)) for c in Fq.c], 't=%.1f' % (time.time()-t0), flush=True)
    pieces = {1: [], 2: [], 3: []}
    for m in range(N+1):
        if m == 0: continue
        pieces[m].append(('ray', Fq.c[m]))
    # endpoint corrections of the integral: sum_{j>=1} deta^j/j! Dtot^{j-1} F  at obs minus at emi
    for sign, dser, where, ev in ((1, dser_o, 'obs', at_obs), (-1, dser_e, 'emi', at_emi)):
        term = Fq; p = Ser.const(1)
        for j in range(1, N+1):
            p = p*dser
            term = Dtot(term, upto=N-j) if j > 1 else trunc(term, N-j)
            add = term*p*sp.Rational(1, sp.factorial(j))
            for m in range(1, N+1):
                if add.c[m] != 0: pieces[m].append((where, sign*ev(add.c[m])))
    if B is not None:
        Bq = pos_expand(B)
        for sign, dser, where, ev in ((1, dser_o, 'obs', at_obs), (-1, dser_e, 'emi', at_emi)):
            Bs = shift_eta(Bq, dser)
            for m in range(1, N+1):
                if Bs.c[m] != 0: pieces[m].append((where, sign*ev(Bs.c[m])))
    return pieces

S_sigma = assemble(I_sigma)
print('S (shear form) assembled, t=%.1f' % (time.time()-t0))
# direct route: ln(Ehat_o/Ehat_e) - ln(That_o/That_e) = -int Gam0kk + [ln(-u.k) - ln That]_e^o  (= -S mod conservation)
S_direct_neg = assemble(I_direct, B_direct)     # this is Lambda^direct = -S^direct
print('Lambda (direct form) assembled, t=%.1f' % (time.time()-t0))

# ---------- classification into named terms ----------
functional_syms = {'k1': set(kfun[1]) | {f_.subs(eta, eta_s) for f_ in kfun[1]},
                   'k2': set(kfun[2]) | {f_.subs(eta, eta_s) for f_ in kfun[2]},
                   'Dx1': set(Dxfun[1]) | {f_.subs(eta, eta_s) for f_ in Dxfun[1]},
                   'Dx2': set(Dxfun[2]) | {f_.subs(eta, eta_s) for f_ in Dxfun[2]},
                   'deta_o1': {deta_o1}, 'deta_o2': {deta_o2}, 'dxo1': set(dxo1), 'dxo2': set(dxo2),
                   'deta_e1': {deta_e1}, 'deta_e2': {deta_e2}}
def field_orders(mono):
    """orders of the field factors of a monomial, with multiplicity (derivatives count as the field)."""
    ords = []
    for fac_, pw in mono.as_powers_dict().items():
        base = fac_.expr if isinstance(fac_, sp.Derivative) else fac_
        if not isinstance(base, sp.core.function.AppliedUndef): continue
        s_ = str(base.func)
        if s_[0] in 'kD': continue
        ords += [int(s_[-1])]*int(pw)
    return ords
def classify(expr):
    groups = {}
    for mono in sp.Add.make_args(sp.expand(expr)):
        atoms = mono.atoms(sp.Symbol) | mono.atoms(sp.core.function.AppliedUndef)
        tags = [kname for kname, ss in functional_syms.items() if atoms & ss]
        # multiplicity of functionals (e.g. Dx1*Dx1)
        mult = []
        for kname in tags:
            deg = sum(sp.degree(mono, f_) for f_ in functional_syms[kname] if mono.has(f_))
            mult.append(kname + ('^%d' % deg if deg > 1 else ''))
        fo = field_orders(mono)
        # how the field functions enter: count factors by order (ignoring derivatives)
        fo_desc = 'fields:' + '+'.join(str(o) for o in sorted(fo, reverse=True)) if fo else 'nofields'
        key = (' * '.join(sorted(mult)) if mult else '-') + ' | ' + fo_desc
        groups.setdefault(key, []).append(mono)
    return {k_: sp.Add(*v_) for k_, v_ in groups.items()}

PHYS = {'k1': 'direction perturbation k1 (aberration + deflection propagated)', 'k2': 'second-order direction k2',
        'Dx1': 'path displacement Dx1 (lensing [transverse] + time delay [radial])', 'Dx2': 'second-order path displacement',
        'deta_o1': 'observer clock shift', 'deta_o2': 'observer clock shift (2nd)', 'dxo1': 'observer displacement',
        'dxo2': 'observer displacement (2nd)', 'deta_e1': 'emission-surface shift', 'deta_e2': 'emission-surface shift (2nd)'}
def term_list(pieces, order, sign=1):
    out = []
    for where, expr in pieces[order]:
        for key, val in classify(sign*expr).items():
            tags = key.split(' | ')[0]
            name = ('[%s] ' % where) + ('; '.join(PHYS.get(t.split('^')[0], t) for t in tags.split(' * ')) if tags != '-' else 'local') + ' | ' + key.split(' | ')[1]
            out.append(dict(name=name, where=where, tags=key, expr=val, nterms=len(sp.Add.make_args(val))))
    return out

L = {}
for m in (1, 2, 3):
    L[m] = term_list(S_sigma, m, sign=-1)       # coefficient of eps^m in Lambda = -S, i.e. L_m/m!
Ld = {}
for m in (1, 2, 3):
    Ld[m] = term_list(S_direct_neg, m, sign=+1)
for m in (1, 2, 3):
    print('=== Lambda coefficient of eps^%d (= L%d/%d!) shear form: %d terms in %d groups' % (m, m, m, sum(t['nterms'] for t in L[m]), len(L[m])))
    for t_ in L[m]: print('   %4d  %s' % (t_['nterms'], t_['name']))
for m in (1, 2, 3):
    print('=== direct form eps^%d: %d terms in %d groups' % (m, sum(t['nterms'] for t in Ld[m]), len(Ld[m])))
    for t_ in Ld[m]: print('   %4d  %s' % (t_['nterms'], t_['name']))

# ---------- m_n and the sky-monopole subtraction (SPEC §7) ----------
# ---------- 'split' display: n = z-hat, k1_3 -> -(Phi1+Psi1) (null condition), Dx1_3 radial ----------
zsub = {n[0]: 0, n[1]: 0, n[2]: 1}
def split_z(expr):
    e = sp.expand(sp.expand(expr).xreplace(zsub))
    e = e.xreplace({kfun[1][2]: -(fname('Phi', 1) + fname('Psi', 1)),
                    kfun[1][2].subs(eta, eta_s): -(fname('Phi', 1) + fname('Psi', 1))})
    return sp.expand(e)
SPLIT_TAGS = {'lens': set([Dxfun[1][0], Dxfun[1][1], Dxfun[1][0].subs(eta, eta_s), Dxfun[1][1].subs(eta, eta_s)]),
              'tdelay': set([Dxfun[1][2], Dxfun[1][2].subs(eta, eta_s)]),
              'k1perp': set([kfun[1][0], kfun[1][1], kfun[1][0].subs(eta, eta_s), kfun[1][1].subs(eta, eta_s)]),
              'lens2': set([Dxfun[2][0], Dxfun[2][1], Dxfun[2][0].subs(eta, eta_s), Dxfun[2][1].subs(eta, eta_s)]),
              'tdelay2': set([Dxfun[2][2], Dxfun[2][2].subs(eta, eta_s)]),
              'k2perp': set([kfun[2][0], kfun[2][1], kfun[2][0].subs(eta, eta_s), kfun[2][1].subs(eta, eta_s)]),
              'k2rad': set([kfun[2][2], kfun[2][2].subs(eta, eta_s)]),
              'deta_o1': {deta_o1}, 'deta_o2': {deta_o2}, 'dxo1': set(dxo1), 'dxo2': set(dxo2), 'deta_e1': {deta_e1}, 'deta_e2': {deta_e2}}
SPLIT_PHYS = {'lens': 'lensing (transverse Dx1)', 'tdelay': 'time delay (radial Dx1)', 'k1perp': 'transverse direction k1 (aberration+deflection)',
              'lens2': 'lensing 2nd (transverse Dx2)', 'tdelay2': 'time delay 2nd (radial Dx2)', 'k2perp': 'transverse k2', 'k2rad': 'radial k2',
              'deta_o1': 'observer clock', 'deta_o2': 'observer clock 2nd', 'dxo1': 'observer displacement', 'dxo2': 'observer displacement 2nd',
              'deta_e1': 'emission-surface shift', 'deta_e2': 'emission-surface shift 2nd'}
def classify_split(expr):
    groups = {}
    for mono in sp.Add.make_args(sp.expand(expr)):
        atoms = mono.atoms(sp.Symbol) | mono.atoms(sp.core.function.AppliedUndef)
        mult = []
        for kname, ss in SPLIT_TAGS.items():
            if atoms & ss:
                deg = sum(sp.degree(mono, f_) for f_ in ss if mono.has(f_))
                mult.append(kname + ('^%d' % deg if deg > 1 else ''))
        fo = field_orders(mono)
        fo_desc = 'fields:' + '+'.join(str(o) for o in sorted(fo, reverse=True)) if fo else 'nofields'
        key = (' * '.join(sorted(mult)) if mult else '-') + ' | ' + fo_desc
        groups.setdefault(key, []).append(mono)
    return {k_: sp.Add(*v_) for k_, v_ in groups.items()}
def term_list_split(pieces, order, sign=1):
    out = []
    for where, expr in pieces[order]:
        for key, val in classify_split(split_z(sign*expr)).items():
            tags = key.split(' | ')[0]
            name = ('[%s] ' % where) + ('; '.join(SPLIT_PHYS.get(t.split('^')[0], t) + ('^'+t.split('^')[1] if '^' in t else '') for t in tags.split(' * ')) if tags != '-' else 'local') + ' | ' + key.split(' | ')[1]
            out.append(dict(name=name, where=where, tags=key, expr=val, nterms=len(sp.Add.make_args(val))))
    return out
Lsplit = {m: term_list_split(S_sigma, m, sign=-1) for m in (1, 2, 3)}
Ldsplit = {m: term_list_split(S_direct_neg, m, sign=+1) for m in (1, 2, 3)}
for m in (1, 2, 3):
    print('=== SPLIT (n=z) shear form eps^%d: %d terms in %d groups' % (m, sum(t['nterms'] for t in Lsplit[m]), len(Lsplit[m])))
    for t_ in sorted(Lsplit[m], key=lambda t: t['name']): print('   %4d  %s' % (t_['nterms'], t_['name']))
for m in (1, 2):
    print('=== SPLIT (n=z) direct form eps^%d: %d terms in %d groups' % (m, sum(t['nterms'] for t in Ldsplit[m]), len(Ldsplit[m])))
    for t_ in sorted(Ldsplit[m], key=lambda t: t['name']): print('   %4d  %s' % (t_['nterms'], t_['name']))
print('--- explicit L2/2 (shear form, n=z) ---')
for t_ in sorted(Lsplit[2], key=lambda t: t['name']): print(t_['name'], ':\n      ', t_['expr'])
print('--- explicit direct-form eps^2 (n=z) ---')
for t_ in sorted(Ldsplit[2], key=lambda t: t['name']): print(t_['name'], ':\n      ', t_['expr'])

L1s, L2s, L3s = sp.symbols('L1 L2 L3')
Lam = Ser([0, L1s, L2s/2, L3s/6])
M = Lam.exp()
m1, m2, m3 = sp.expand(M.c[1]), sp.expand(2*M.c[2]), sp.expand(6*M.c[3])
a1, a2, a3 = sp.symbols('<m1> <m2> <m3>')
Mbar = Ser([1, a1, a2/2, a3/6])
Th = (M/Mbar)
Theta = {1: sp.expand(Th.c[1]), 2: sp.expand(2*Th.c[2]), 3: sp.expand(6*Th.c[3])}
print('m1 =', m1, '; m2 =', m2, '; m3 =', m3)
print('Theta1 =', Theta[1]); print('Theta2 =', Theta[2]); print('Theta3 =', Theta[3])

with open('L_terms.pkl', 'wb') as f:
    pickle.dump(dict(doc=__doc__, L=L, L_direct=Ld, L_split=Lsplit, L_direct_split=Ldsplit, emission=emission, definitions=ray.definitions,
                     S_pieces=S_sigma, Lambda_direct_pieces=S_direct_neg, m=(m1, m2, m3), Theta=Theta,
                     symbols=dict(n=n, eta_star=eta_s, eta_o_bar=eta_o, deta_o=(deta_o1, deta_o2), dxo=(dxo1, dxo2),
                                  deta_e=(deta_e1, deta_e2), kfun=kfun, Dxfun=Dxfun)), f)
print('saved L_terms.pkl, t=%.1f' % (time.time()-t0))
