"""
isw_structure.py -- agent C5: general-l STRUCTURE of the closed-form ISW integral
   I_l(k) = (4/sqrt3) int_{s_*}^{S} [j_2(s/sqrt3)/(s/sqrt3)] j_l(S-s) ds ,  s_* = k eta_*, S = k eta_o, w := S - s_*,
and its explicit forms for l = 0..6 (written to isw_closed_forms.tex).

THEOREM (general l; derivation in report_C5.md / linear_cl.tex):
   I_l = E_l(s_*, S) + sum_{sigma=±} { A^sigma_l(S) [Ci(a_sigma S) - Ci(a_sigma s_*)] + B^sigma_l(S) [Si(a_sigma S) - Si(a_sigma s_*)]
                                      + C^sigma_l(S) Cin(a_sigma w)  + D^sigma_l(S) Si(a_sigma w) },
   a_± = 1 ± 1/sqrt3,  Cin(z) = int_0^z (1-cos t)/t dt = EulerGamma + ln z - Ci(z) (entire),
   E_l elementary (Laurent polynomials in s_*, S, w times sin/cos of s_*/sqrt3, S - s_*, S/sqrt3, S).
 With h = h_l^(1) (spherical Hankel), g(x) = h_2^(1)(x)/x, D = d/dS:
   A^+ - i B^+ = (2/sqrt3) conj[ P(D - i) h_l(S) ] ... see residue formulas below (all derived from
   integrand = (1/sqrt3) (g + gbar)(x = s/sqrt3) (h + hbar)(S - s) and  int e^{i b s} ds/s = Ci(|b|s) + i sgn(b) Si(|b|s)):
   s-endpoint (pole of g's Laurent factor at s = 0, order <= 4  ->  derivatives of h_l up to 3rd order):
      rho_{gH}    = Res_{s=0} Lg(s/sqrt3) LH(S-s),     A^- = (2/sqrt3) Re[e^{iS} rho_gH],    B^- = (2/sqrt3) Im[e^{iS} rho_gH]
      rho_{gHbar} = Res_{s=0} Lg(s/sqrt3) LHbar(S-s),  A^+ = (2/sqrt3) Re[e^{-iS} rho_gHb],  B^+ = -(2/sqrt3) Im[e^{-iS} rho_gHb]
      equivalently e^{iS} rho_gH = i[-3(D-i) + (9 sqrt3/2) i (D-i)^2 + (9/2)(D-i)^3] h_l^(1)(S)  (and c.c. with D+i, h^(2)).
   u-endpoint (pole of h_l's Laurent factor at u = S - s = 0, order <= l+1  ->  derivatives of g up to order l):
      tau_{gH}    = Res_{u=0} LH(u) Lg((S-u)/sqrt3),    C^- = -(2/sqrt3) Re[e^{iS/sqrt3} tau_gH],  D^- = -(2/sqrt3) Im[e^{iS/sqrt3} tau_gH]
      tau_{gHbar} = Res_{u=0} LHbar(u) Lg((S-u)/sqrt3), C^+ = -(2/sqrt3) Re[e^{iS/sqrt3} tau_gHb], D^+ = (2/sqrt3) Im[e^{iS/sqrt3} tau_gHb]
   where h_l^(1)(u) = e^{iu} LH(u), LH(u) = (-i)^{l+1} sum_{k=0}^{l} i^k (l+k)!/(k!(l-k)! 2^k) u^{-k-1}, LHbar = conj coefficients,
   g(x) = e^{ix} Lg(x), Lg(x) = i (x^-2 + 3 i x^-3 - 3 x^-4).
 Regularity at the observer (u -> 0) <=> Re[e^{iS/sqrt3}(tau_gH + tau_gHb)] = 0 (the ln w coefficient): checked.
 Observed symmetry (checked l <= 6): A^- = -A^+, B^- = -B^+, C^- = -C^+, D^- = D^+.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
import mpmath as mp
from linear_cl import I_l_closed, I_l_closed_eta0, S, sst, A_PLUS, A_MINUS, SQ3, HERE

I = sp.I
v = sp.Symbol('v')          # local variable (s or u)
w = sp.Symbol('w', positive=True)
TR = sp.Symbol('t', real=True)


def LH(l, uu, conj=False):
    c = sum(I**kk * sp.factorial(l + kk) / (sp.factorial(kk) * sp.factorial(l - kk) * 2**kk) * uu**(-kk - 1) for kk in range(l + 1))
    e = (-I)**(l + 1) * c
    return sp.conjugate(e) if conj else e


def Lg(xx):
    return I * (xx**-2 + 3 * I * xx**-3 - 3 * xx**-4)


def _res(expr, var):
    return sp.expand(sp.series(expr, var, 0, 1).removeO().coeff(var, -1)) if True else None


def _residue(expr, var):
    return sp.residue(sp.expand(expr), var, 0)


def coefficients(l):
    """Predicted A^±, B^±, C^±, D^± from the residue formulas (sympy, real expressions in S).
       The Ci/Si coefficient of a phase term e^{i b s} R(s) is Res_{s=0}[e^{i b s} R(s)] (higher poles feed Ci/Si through
       the recursion int e^{ibs}/s^m = ... + (ib)^{m-1}/(m-1)! int e^{ibs}/s), hence residues of the FULL functions:
         e^{iS} rho_gH     = Res_{s=0} g(s/sqrt3) h^(1)_l(S-s),   e^{-iS} rho_gHb = Res_{s=0} g(s/sqrt3) h^(2)_l(S-s),
         e^{iS/sqrt3} tau_gH = Res_{u=0} h^(1)_l(u) g((S-u)/sqrt3), e^{iS/sqrt3} tau_gHb = Res_{u=0} h^(2)_l(u) g((S-u)/sqrt3)."""
    Sr = sp.Symbol('S', real=True)
    # residue of (Laurent polynomial L(v)) x (function A(v) analytic at 0) = sum_j [v^-j]L * A^(j-1)(0)/(j-1)!
    # (sympy's residue() was found to be unreliable here, so it is done by hand)
    def res(L, A):
        out = 0
        for t in sp.Add.make_args(sp.expand(L)):
            c, e = t.as_coeff_exponent(v)
            j = -int(e)
            if j >= 1:
                out += c * sp.diff(A, v, j - 1).subs(v, 0) / sp.factorial(j - 1)
        return sp.expand(out)
    LHb = sp.conjugate(LH(l, TR)).subs(TR, v)
    Lgv = sp.expand(Lg(v / SQ3))
    rho1 = res(Lgv, sp.exp(I * v / SQ3) * sp.exp(I * (Sr - v)) * LH(l, Sr - v))                      # Res g(s/sqrt3) h1(S-s)
    rho2 = res(Lgv, sp.exp(I * v / SQ3) * sp.exp(-I * (Sr - v)) * LHb.subs(v, Sr - v))              # Res g(s/sqrt3) h2(S-s)
    gS = sp.exp(I * (Sr - v) / SQ3) * Lg((Sr - v) / SQ3)
    tau1 = res(LH(l, v), sp.exp(I * v) * gS)                                                          # Res h1(u) g((S-u)/sqrt3)
    tau2 = res(LHb, sp.exp(-I * v) * gS)                                                              # Res h2(u) g((S-u)/sqrt3)
    two = 2 / SQ3
    re = lambda z: sp.simplify(sp.re(sp.expand_complex(sp.expand(z))))
    im = lambda z: sp.simplify(sp.im(sp.expand_complex(sp.expand(z))))
    out = {
        'A-': two * re(rho1), 'B-': two * im(rho1),
        'A+': two * re(rho2), 'B+': -two * im(rho2),
        'C-': -two * re(tau1), 'D-': -two * im(tau1),
        'C+': -two * re(tau2), 'D+': two * im(tau2),
        'logw': re(tau1 + tau2),
    }
    return {kk: sp.simplify(vv.subs(Sr, S)) for kk, vv in out.items()}


def real_forms(l):
    """General-l real forms of the transcendental coefficients ('+' members; '-' members by the symmetries):
         A_l(S) = (2/sqrt3) Res_{s=0} Y(s) y_l(S-s) = sqrt3 [3 y_l^(3)(S) + y_l'(S)]
         B_l(S) = -(2/sqrt3) Res_{s=0} Y(s) j_l(S-s) = -sqrt3 [3 j_l^(3)(S) + j_l'(S)]
         C_l(S) = -(2/sqrt3) Res_{u=0} y_l(u) Y(S-u),   D_l(S) = -(2/sqrt3) Res_{u=0} y_l(u) J(S-u),
       Y(s) = y_2(s/sqrt3)/(s/sqrt3) = -27/s^4 - 3/(2 s^2) + regular,  J(s) = j_2(s/sqrt3)/(s/sqrt3) (entire)."""
    z = sp.Symbol('z')
    yl = sp.expand_func(sp.yn(l, z)); jl = sp.expand_func(sp.jn(l, z))
    A = SQ3 * (3 * sp.diff(yl, z, 3) + sp.diff(yl, z)).subs(z, S)
    B = -SQ3 * (3 * sp.diff(jl, z, 3) + sp.diff(jl, z)).subs(z, S)
    Yf = (sp.expand_func(sp.yn(2, z)) / z).subs(z, (S - v) / SQ3)
    Jf = (sp.expand_func(sp.jn(2, z)) / z).subs(z, (S - v) / SQ3)
    pp = sp.series(sp.expand_func(sp.yn(l, v)), v, 0, 1).removeO()     # principal part of y_l (pole order l+1)
    def res(L, F):
        out = 0
        for t in sp.Add.make_args(sp.expand(L)):
            c, e = t.as_coeff_exponent(v)
            j = -int(e)
            if j >= 1:
                out += c * sp.diff(F, v, j - 1).subs(v, 0) / sp.factorial(j - 1)
        return out
    C = -2 / SQ3 * res(pp, Yf)
    D = -2 / SQ3 * res(pp, Jf)
    return {'A': sp.simplify(A), 'B': sp.simplify(B), 'C': sp.simplify(C), 'D': sp.simplify(D)}


def _match(atom):
    """Classify a Ci/Si atom: returns (kind, sigma, endpoint) with endpoint in {'S','s*','w'}."""
    arg = atom.args[0]
    for sig, a in (('+', A_PLUS), ('-', A_MINUS)):
        for ep, val in (('S', S), ('s*', sst), ('w', S - sst)):
            if sp.simplify(sp.expand(arg - a * val)) == 0:
                return (type(atom).__name__, sig, ep)
    raise ValueError(atom)


def structure(l):
    """Transcendental part T_l built from the predicted coefficients, and E_l = I_l - T_l."""
    co = coefficients(l)
    T = 0
    for sig, a in (('+', A_PLUS), ('-', A_MINUS)):
        T += co['A' + sig] * (sp.Ci(a * S) - sp.Ci(a * sst)) + co['B' + sig] * (sp.Si(a * S) - sp.Si(a * sst))
        Cin = sp.EulerGamma + sp.log(a) + sp.log(S - sst) - sp.Ci(a * (S - sst))
        T += co['C' + sig] * Cin + co['D' + sig] * sp.Si(a * (S - sst))
    return co, T


def verify(l, log=print):
    t0 = time.time()
    e = I_l_closed(l)
    co, T = structure(l)
    ok = True
    # (1) coefficient-by-coefficient comparison of every Ci/Si atom of the brute-force closed form (numerical, 40 digits)
    pts = [(sp.Rational(3, 10), sp.Rational(27, 10)), (sp.Rational(7, 5), sp.Rational(53, 7)), (sp.Rational(1, 7), sp.Rational(9, 2))]
    worst = 0
    for atom in e.atoms(sp.Ci, sp.Si):
        kind, sig, ep = _match(atom)
        c_brute = e.coeff(atom)
        if ep == 'S':
            pred = co[('A' if kind == 'Ci' else 'B') + sig]
        elif ep == 's*':
            pred = -co[('A' if kind == 'Ci' else 'B') + sig]
        else:
            pred = -co['C' + sig] if kind == 'Ci' else co['D' + sig]
        for (a_, b_) in pts:
            d = sp.N((c_brute - pred).subs({sst: a_, S: b_}), 40)
            worst = max(worst, abs(d))
    ok &= worst < 1e-30
    # (2) regularity: the ln w coefficient vanishes identically
    lw = sp.simplify(co['logw'])
    ok &= (lw == 0)
    # (3) symmetries
    sym = all(sp.simplify(co[a] + co[b]) == 0 for a, b in (('A-', 'A+'), ('B-', 'B+'), ('C-', 'C+'))) and sp.simplify(co['D-'] - co['D+']) == 0
    # (3b) real residue forms: A = sqrt3 (3 y_l''' + y_l')(S), B = -sqrt3 (3 j_l''' + j_l')(S),
    #      C = -(2/sqrt3) Res_{u=0} y_l(u) Y(S-u), D = -(2/sqrt3) Res_{u=0} y_l(u) J(S-u)
    rf = real_forms(l)
    for key, kk in (('A', 'A+'), ('B', 'B+'), ('C', 'C+'), ('D', 'D+')):
        for (a_, b_) in pts:
            worst = max(worst, abs(sp.N((rf[key] - co[kk]).subs(S, b_), 40)))
    ok &= worst < 1e-30
    # (4) E_l = I_l - T elementary: no Ci/Si after expansion, and logs/EulerGamma cancel numerically
    E = sp.expand(e - sp.expand(T))
    left = E.atoms(sp.Ci, sp.Si)
    ok &= (len(left) == 0)
    Ef = E.subs({sp.EulerGamma: 0})
    E_noLog = Ef.replace(lambda z: isinstance(z, sp.log), lambda z: sp.Integer(0))
    logpart = 0
    for (a_, b_) in pts:
        logpart = max(logpart, abs(sp.N((E - E_noLog).subs({sst: a_, S: b_}), 40)))
    ok &= logpart < 1e-30
    log("  l=%d: Ci/Si coefficients (residue formulas) vs brute force: max diff %.1e; ln(w) coefficient = %s; "
        "symmetries A-=-A+, B-=-B+, C-=-C+, D-=D+: %s; real residue forms incl. in max diff; E_l elementary (no Ci/Si, logs cancel %.0e) [%.0fs]"
        % (l, worst, lw, sym, logpart, time.time() - t0))
    return ok, co, sp.expand(E_noLog)


def E_collected(E):
    """Collect E_l by its trigonometric factors."""
    trig = sorted({f for t in sp.Add.make_args(E) for f in sp.Mul.make_args(t) if isinstance(f, (sp.sin, sp.cos))
                   or (f.is_Pow and isinstance(f.base, (sp.sin, sp.cos)))}, key=str)
    return sp.collect(sp.expand(sp.expand_trig(E) if False else E), trig, evaluate=False)


def write_tex(lmax=6, log=print):
    lines = [r"% auto-generated by isw_structure.py (agent C5).  $s_*=k\eta_*$, $S=k\eta_o$, $w=S-s_*$, $a_\pm=1\pm1/\sqrt3$.",
             r"% $I_\ell=E_\ell+A_\ell\,\mathcal{C}+B_\ell\,\mathcal{S}+C_\ell[\mathrm{Cin}(a_+w)-\mathrm{Cin}(a_-w)]+D_\ell[\mathrm{Si}(a_+w)+\mathrm{Si}(a_-w)]$,",
             r"% $\mathcal{C}=\mathrm{Ci}(a_+S)-\mathrm{Ci}(a_+s_*)-\mathrm{Ci}(a_-S)+\mathrm{Ci}(a_-s_*)$, $\mathcal{S}$ likewise with Si."]
    allok = True
    for l in range(0, lmax + 1):
        ok, co, E = verify(l, log)
        allok &= ok
        rf = real_forms(l)
        lines.append(r"\paragraph{$\ell=%d$.}" % l)
        for key in ('A', 'B', 'C', 'D'):
            lines.append(r"\begin{dmath*} %s_{%d}(S) = %s \end{dmath*}" % (key, l, sp.latex(sp.factor(rf[key]))))
        coll = E_collected(E)
        parts = []
        for trig, c in coll.items():
            parts.append(r"\left(%s\right) %s" % (sp.latex(sp.factor(sp.together(c))), sp.latex(trig)))
        lines.append(r"\begin{dmath*} E_{%d} = %s \end{dmath*}" % (l, " + ".join(parts)))
        log("     E_%d: %d terms, %d trigonometric structures" % (l, len(sp.Add.make_args(E)), len(coll)))
    with open(os.path.join(HERE, 'isw_closed_forms.tex'), 'w') as fh:
        fh.write("\n".join(lines) + "\n")
    return allok


if __name__ == '__main__':
    ok = write_tex()
    print("STRUCTURE ALL PASS:", ok)


def write_coeff_tex(lmax=6):
    """Short file: only the transcendental coefficients A_l..D_l (l = 0..lmax), for the LaTeX section."""
    lines = [r"% auto-generated by isw_structure.write_coeff_tex (agent C5); $S=k\eta_o$, $z=S/\sqrt3$."]
    z = sp.Symbol('z', positive=True)
    for l in range(lmax + 1):
        rf = real_forms(l)
        lines.append(r"\paragraph{$\ell=%d$.}" % l)
        for key in ('A', 'B', 'C', 'D'):
            e = sp.factor(rf[key])
            if key in 'CD':
                e = sp.factor(sp.simplify(e.subs(S, SQ3 * z)))
            lines.append(r"\begin{dmath*} %s_{%d} = %s \end{dmath*}" % (key, l, sp.latex(e)))
    with open(os.path.join(HERE, 'isw_coeffs.tex'), 'w') as fh:
        fh.write("\n".join(lines) + "\n")
