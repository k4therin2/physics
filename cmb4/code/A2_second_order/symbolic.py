"""symbolic.py -- explicit K_2(k1,k2;n) channel kernels in NAMED OBJECTS (sympy), generated from C2's direct-form L2 list.

Fourier rules (C2 convention, x_o = 0, ray xbar = n chi, chi = eta_o - eta; f(x) = int d^3k/(2pi)^3 f_k e^{ik.x}):
  a first-order field X1 (zeta_k1) on the ray -> T_X(k1, eta) E1,  E1 := e^{i k1 mu1 chi},  mu1 = n.khat1 (likewise slot 2);
  d_eta -> eta-derivative of the transfer function (symbol suffix 'p', 'pp'), d_j -> i k_j;  v1^j = i k^j T_V.
  Geometry: transverse products k1perp.k2perp = k1 k2 (mu - mu1 mu2), mu = khat1.khat2, K = k1 + k2, Kn := n.K = k1 mu1 + k2 mu2,
  K.k1 = k1^2 + k1 k2 mu, K.k2 = k2^2 + k1 k2 mu.
NAMED FIRST-ORDER OBJECTS (symbols carry slot index 1,2 and are evaluated at the term's time: eta on the ray, eta_o at 'obs',
eta_* at 'emi'):  TPhi, TPhip (= d T_Phi/d eta), TdT, TdTp, TV, TVp;  ray functionals (C2 'definitions' with the corrected sign):
  Kf_a(eta)  := T_V(k_a,eta_o) + 2 int_eta^{eta_o} T_Phi(k_a,eta') E_a(eta') deta'        [k1perp = i k_perp Kf  : aberration + deflection]
  Dp_a(eta)  := X_o(k_a) - int_eta^{eta_o} Kf_a(eta') deta'                                [Dx1perp = i k_perp Dp : lensing displacement]
  Dz_a(eta)  := i k_a mu_a X_o(k_a) + deta_o1(k_a) + 2 int_eta^{eta_o} T_Phi E_a deta'      [n.Dx1 : time delay]
  X_o(k)     := int_0^{eta_o} T_V(k,eta) deta  (observer displacement dx_o1 = i k X_o),  deta_o1(k) := -(1/eta_o) int_0^{eta_o} eta T_Phi deta  (clock),
  deta_e1    := eta_* T_dT(k,eta_*) E(eta_*)   (emission-surface shift).
SECOND-ORDER OBJECTS (pair kernels of C1, functions of (k1,k2,K,eta); phase EK := e^{i Kn chi} = E1 E2):
  Phi2, Phi2p, Psi2, Psi2p, dT2, V2 (v2 = grad V2), KvT12 := K_vT(k1,k2), KvT21, Kom12, Kom21, H, Hp (h2_ij = Lambda_ij,kl(K) k1_k k1_l H).
  n.v2T = (1/2)[n.P(K)k1 KvT12 + n.P(K)k2 KvT21],  n.P(K)k_a = k_a mu_a - Kn (K.k_a)/K^2;  om likewise;
  h_nn = [(n.P k1)^2 - (1/2)(1 - Kn^2/K^2)(k1^2 - (K.k1)^2/K^2)] H;  (h_11 + h_22) = -h_nn (traceless).
K_2 = 2 x [symmetrised kernel of C2's L2/2 list] + K1(k1)K1(k2) - Delta_0(k1)K1(k2) - Delta_0(k2)K1(k1).
Outputs: symbolic_kernels.pkl (channel -> where, expr) and symbolic_kernels.tex.  __main__ verifies numerically against kernel2."""
import os, sys, re, pickle
import numpy as np
import sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__))
LT = pickle.load(open(os.path.join(HERE, '..', 'C2_observable', 'L_terms.pkl'), 'rb'))
ETA, XS, YS, ZS = sp.symbols('eta x y z')
k1, k2, mu, mu1, mu2, K, Kn, etas, etao = sp.symbols('k1 k2 mu mu1 mu2 K Kn eta_star eta_o', real=True)
Pperp = sp.Symbol('Pperp')          # k1perp.k2perp = k1 k2 (mu - mu1 mu2)
E = {1: sp.Symbol('E1'), 2: sp.Symbol('E2')}
a = {1: sp.symbols('a1 a2 a3'), 2: sp.symbols('b1 b2 b3')}
kmag = {1: k1, 2: k2}; kmu = {1: mu1, 2: mu2}
def S(name, s): return sp.Symbol('%s_%d' % (name, s))
second = {nm: sp.Symbol(nm) for nm in ['Phi2', 'Phi2p', 'Psi2', 'Psi2p', 'dT2', 'V2', 'KvT12', 'KvT21', 'Kom12', 'Kom21', 'H', 'Hp']}
EK = sp.Symbol('EK')
_re_scalar = re.compile(r'^(Psi|Phi|dT)(\d)$'); _re_vec = re.compile(r'^(v|om)_(\d)_(\d)$'); _re_h = re.compile(r'^h_(\d)(\d)_(\d)$')
_re_fun = re.compile(r'^(k1|Dx1)_(\d)$')
TNAME = {'Psi': 'TPhi', 'Phi': 'TPhi', 'dT': 'TdT'}

def first_obj(f, s):
    """sympy image of a first-order factor for slot s (phase included)."""
    der = [0, 0, 0, 0]
    if isinstance(f, sp.Derivative):
        for v_, c_ in f.variable_count: der[(ETA, XS, YS, ZS).index(v_)] += int(c_)
        f = f.expr
    if isinstance(f, sp.Symbol):
        nm = f.name
        if nm == 'deta_o1': return S('detao1', s)
        if nm == 'deta_e1': return S('detae1', s)
        if nm.startswith('dxo1_'):
            c = int(nm[-1]) - 1
            return sp.I * a[s][c] * S('Xo', s)
        raise ValueError(nm)
    nm = f.func.__name__
    m = _re_fun.match(nm)
    if m:
        c = int(m.group(2)) - 1
        if m.group(1) == 'k1':
            return sp.I * a[s][c] * S('Kf', s) if c < 2 else -2 * S('TPhi', s) * E[s]
        return sp.I * a[s][c] * S('Dp', s) if c < 2 else S('Dz', s)
    fac = (sp.I * a[s][0])**der[1] * (sp.I * a[s][1])**der[2] * (sp.I * a[s][2])**der[3] * E[s]
    m = _re_scalar.match(nm)
    if m:
        return S(TNAME[m.group(1)] + 'p' * der[0], s) * fac
    m = _re_vec.match(nm)
    if m and m.group(1) == 'v':
        c = int(m.group(2)) - 1
        return sp.I * a[s][c] * S('TV' + 'p' * der[0], s) * fac
    raise ValueError(nm)

def nPk(s):
    Kk = kmag[s]**2 + k1 * k2 * mu
    return kmag[s] * kmu[s] - Kn * Kk / K**2

def second_obj(f):
    der = [0, 0, 0, 0]
    if isinstance(f, sp.Derivative):
        for v_, c_ in f.variable_count: der[(ETA, XS, YS, ZS).index(v_)] += int(c_)
        f = f.expr
    nm = f.func.__name__
    assert der[1] == der[2] == 0, f
    fac = (sp.I * Kn)**der[3] * EK
    m = _re_scalar.match(nm)
    if m:
        base = {'Phi': 'Phi2', 'Psi': 'Psi2', 'dT': 'dT2'}[m.group(1)] + 'p' * der[0]
        return second[base] * fac
    m = _re_vec.match(nm)
    if m:
        assert int(m.group(2)) == 3 and der[0] == 0, f
        if m.group(1) == 'v':
            return (sp.I * Kn * second['V2'] + (nPk(1) * second['KvT12'] + nPk(2) * second['KvT21']) / 2) * fac
        return ((nPk(1) * second['Kom12'] + nPk(2) * second['Kom21']) / 2) * fac
    m = _re_h.match(nm)
    if m:
        i, j = int(m.group(1)), int(m.group(2))
        assert i == j and i in (1, 2), f
        hnn = (nPk(1)**2 - sp.Rational(1, 2) * (1 - Kn**2 / K**2) * (k1**2 - (k1**2 + k1 * k2 * mu)**2 / K**2)) * second['H' + 'p' * der[0]]
        return -hnn / 2 * fac        # h_11 + h_22 = -h_nn ; each of h_11, h_22 contributes half of the pair
    raise ValueError(nm)

def is_second(f):
    g = f.expr if isinstance(f, sp.Derivative) else f
    return isinstance(g, sp.core.function.AppliedUndef) and g.func.__name__[-1] == '2' and not _re_fun.match(g.func.__name__)

def reduce_perp(e):
    e = sp.expand(e.subs({a[1][1]: 0, a[2][1]: 0}))
    e = e.subs(a[1][0] * a[2][0], Pperp)
    assert not e.has(a[1][0]) and not e.has(a[2][0]), e
    return e.subs({a[1][2]: k1 * mu1, a[2][2]: k2 * mu2})

def channel_kernels(form='L_direct_split'):
    out = {}
    for term in LT[form][2]:
        where = term['where']; name = term['name']
        tot = 0
        for mono in sp.Add.make_args(sp.expand(term['expr'])):
            coeff = 1; facs = []
            for f in sp.Mul.make_args(mono):
                base, ex = (f.base, int(f.exp)) if isinstance(f, sp.Pow) else (f, 1)
                if base.is_number: coeff *= base**ex; continue
                if base == ETA or base == sp.Symbol('eta'): coeff *= ETA**ex; continue
                if base == sp.Symbol('eta_star'): coeff *= etas**ex; continue
                facs += [base] * ex
            if len(facs) == 1 and is_second(facs[0]):
                if facs[0].has(sp.Function('h_11_2')) or facs[0].has(sp.Function('h_22_2')):
                    tot += coeff * second_obj(facs[0])      # each of h11, h22 gives -h_nn/2
                else:
                    tot += coeff * second_obj(facs[0])
            else:
                assert len(facs) == 2, facs
                A, B = facs
                tot += coeff * (first_obj(A, 1) * first_obj(B, 2) + first_obj(A, 2) * first_obj(B, 1)) / 2
        tot = 2 * reduce_perp(tot)          # K_2 = 2 x kernel of L2/2
        key = name
        if key in out: out[key] = (where, out[key][1] + tot)
        else: out[key] = (where, tot)
    return {k_: (w, sp.simplify(sp.factor_terms(e))) for k_, (w, e) in out.items()}

if __name__ == '__main__':
    ch = channel_kernels()
    pickle.dump({'channels': {k_: (w, sp.srepr(e)) for k_, (w, e) in ch.items()}, 'doc': __doc__}, open(os.path.join(HERE, 'symbolic_kernels.pkl'), 'wb'))
    with open(os.path.join(HERE, 'symbolic_kernels.tex'), 'w') as fh:
        fh.write('% K_2 channel kernels (direct form), generated by symbolic.py. ray: int_{eta_*}^{eta_o} d eta of the expression.\n')
        fh.write('% Pperp = k1 k2 (mu - mu1 mu2); E_a = exp(i k_a mu_a chi); EK = E1 E2.\n')
        for k_, (w, e) in ch.items():
            fh.write('\\paragraph{%s [%s]}\n\\begin{dmath*}\n%s\n\\end{dmath*}\n' % (k_.replace('|', '/').replace('_', '\\_').replace('<', '$<$').replace('>', '$>$'), w, sp.latex(e)))
        fh.write('\\paragraph{B channels} $K_1(k_1)K_1(k_2) - \\Delta_0(k_1)K_1(k_2) - \\Delta_0(k_2)K_1(k_1)$.\n')
    for k_, (w, e) in ch.items():
        print('[%s] %s :\n    %s' % (w, k_, e))
