"""X1: internal consistency of Adam's two printed forms of the O(k0^6) superhorizon quadrupole coefficients (his Sec. 9.2):
explicit polynomials vs the three-tensor form with beta_H, beta_a (eq. (85) and the display after it)."""
import sympy as sp
o, e, i, t = sp.symbols('o e i t', positive=True)
D = o**2 - e**2
L = sp.integrate(t * (o - t)**2, (t, e, o)); E = o**4 - e**4
bH = -sp.Rational(11, 14) * L - sp.Rational(3, 40) * E
bH_poly = -(o - e) / 420 * (59*o**3 + 59*o**2*e - 106*o*e**2 + 114*e**3)
ba = sp.Rational(27, 7) * L + sp.Rational(7, 20) * E + (o**4 + 3*i**4) / 4 - e**4
ba_poly = sp.Rational(3, 140) * (43*o**4 - 90*o**2*e**2 + 120*o*e**3 - 108*e**4 + 35*i**4)
print('beta_H forms agree:', sp.simplify(bH - bH_poly) == 0, '| beta_a forms agree:', sp.simplify(ba - ba_poly) == 0)
c11_poly = -4*sp.pi/637875 * (o-e)**2*(o+e)*(59*o**3+59*o**2*e-106*o*e**2+114*e**3)
c11_tens = 16*sp.pi/6075 * D * bH
q22_poly = -8*sp.pi/637875 * D * (105*i**4+719*o**4-1920*o**2*e**2+2560*o*e**3-1464*e**4)
q22_tens = 32*sp.pi/18225 * D * (30*bH - ba)
q13_poly = -2*sp.pi/1913625 * D * (2835*i**4+9973*o**4-25440*o**2*e**2+33920*o*e**3-21288*e**4)
q13_tens = 8*sp.pi/18225 * D * (110*bH - 9*ba)
tot_poly = -2*sp.pi/1913625 * D * (4095*i**4+18601*o**4-48480*o**2*e**2+64640*o*e**3-38856*e**4)
tot_tens = 8*sp.pi/18225 * D * (230*bH - 13*ba)
for nm, a, b in (('C11 k^6', c11_poly, c11_tens), ('1/4C22 k^6', q22_poly, q22_tens), ('1/3ReC13 k^6', q13_poly, q13_tens),
                 ('total eps^4 k^6', tot_poly, tot_tens), ('22+13 (poly) vs total (poly)', q22_poly + q13_poly, tot_poly)):
    print('%-30s poly - tensor form = %s' % (nm, sp.factor(sp.simplify(a - b))))
sub = {o: 8, e: 2, i: 0}
print('at o=8,e=2,i=0: q13 poly %.4f, q13 tensor %.4f, q22 poly %.4f, q22 tensor %.4f' % tuple(float(v.subs(sub)) for v in (q13_poly, q13_tens, q22_poly, q22_tens)))
