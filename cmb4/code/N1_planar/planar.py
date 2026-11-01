"""
N1_planar: exact (nonperturbative) planar-symmetric GR + radiation fluid, comoving-orthogonal gauge,
pseudo-spectral in x (period 2 pi, k = 1), DOP853 in time, exact ray tracing to a comoving observer.

Conventions (see README.md):  k = 1; a(t) = t (background conformal time); That_bar = 1 (conformal temperature);
T_phys = 1/(t alh); zeta(x) = eps cos x;  T_* = 1/t_star;  observer = fluid element x = x_o;
observation at proper time Delta tau = (t_obar^2 - t_i^2)/2 after crossing T_phys = T_i = 1/t_i.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
import gen_numeric as gn


class PlanarRun:
    def __init__(self, eps, N=64, t_i=0.01, t_star=2.0, t_obar=8.0, x_o=0.0, cH=-0.5,
                 rtol=1e-13, atol=1e-16, ic_order=6, nmu=64, verbose=True):
        self.eps, self.N, self.t_i, self.t_star, self.t_obar, self.x_o = eps, N, t_i, t_star, t_obar, x_o
        self.cH, self.rtol, self.atol, self.ic_order, self.nmu, self.verbose = cH, rtol, atol, ic_order, nmu, verbose
        self.x = 2*np.pi*np.arange(N)/N
        self.m = np.arange(N//2 + 1)
        self.ik = 1j*self.m
        self.ik[-1] = 0.0            # Nyquist derivative -> 0
        self.ics = gn.load_ics(ic_order)
        self.log = []

    # ---------- spectral helpers ----------
    def zeta_derivs(self, x, nmax):
        return [self.eps*np.cos(x + n*np.pi/2) for n in range(nmax + 1)]

    def dx(self, f, n=1):
        return np.fft.irfft(np.fft.rfft(f)*self.ik**n, n=self.N)

    def interp_setup(self, xr):
        """matrix for evaluating rfft coefficients at arbitrary points xr (real fields)."""
        E = np.exp(1j*np.outer(xr, self.m))/self.N
        E[:, 1:-1] *= 2.0
        E[:, -1] = np.cos(self.N/2*xr)/self.N    # Nyquist mode: symmetric choice
        return E

    @staticmethod
    def interp(E, fhat):
        return (E @ fhat).real if E.shape[1] == fhat.shape[0] else None

    # ---------- initial data (gradient expansion of the regular solution) ----------
    def series_fields(self, t, x):
        """Ah, Bh, pA, pB, alh at time t and points x from the gradient expansion (nonlinear in zeta)."""
        zs = self.zeta_derivs(x, self.ic_order)
        u = np.zeros_like(x); v = np.zeros_like(x); ut = np.zeros_like(x); vt = np.zeros_like(x)
        for n, (fa, fb) in self.ics.items():
            an = fa(*zs) + 0*x; bn = fb(*zs) + 0*x
            u += an*t**n; v += bn*t**n; ut += n*an*t**(n - 1); vt += n*bn*t**(n - 1)
        Ah = np.exp(zs[0])*(1 + u); Bh = np.exp(zs[0])*(1 + v)
        alh = ((1 + u)*(1 + v)**2)**(1/3)
        pA = ut/(1 + u)/alh; pB = vt/(1 + v)/alh
        return Ah, Bh, pA, pB, alh

    def initial_state(self):
        Ah, Bh, pA, pB, alh = self.series_fields(self.t_i, self.x)
        return np.concatenate([Ah, Bh, pA, pB])

    # ---------- field equations ----------
    def unpack(self, y):
        N = self.N
        return y[:N], y[N:2*N], y[2*N:3*N], y[3*N:4*N]

    def field_quantities(self, t, y):
        Ah, Bh, pA, pB = self.unpack(y)
        z0, z1, z2 = self.zeta_derivs(self.x, 2)
        larg = (np.log(Ah) + 2*np.log(Bh))/3 - z0
        alh = np.exp(larg); dalh = np.expm1(larg)
        Ax, Axx = self.dx(Ah), self.dx(Ah, 2)
        Bx, Bxx = self.dx(Bh), self.dx(Bh, 2)
        pAx, pBx = self.dx(pA), self.dx(pB)
        return (t, Ah, Bh, alh, dalh, Ax, Bx, Axx, Bxx, pA, pB, pAx, pBx, z0, z1, z2)

    def rhs_fields(self, t, y):
        q = self.field_quantities(t, y)
        Ah, Bh, alh, pA, pB = q[1], q[2], q[3], q[9], q[10]
        dAh = alh*Ah*pA
        dBh = alh*Bh*pB
        dpA = gn.f_pA_t(*q)
        dpB = gn.f_pB_t(*q)
        if self.cH != 0.0:
            dpB = dpB + self.cH*t**2*gn.f_Ham(*q)
        return np.concatenate([dAh, dBh, dpA, dpB])

    def constraints(self, t, y):
        q = self.field_quantities(t, y)
        H = gn.f_Ham(*q); M = gn.f_Mom(*q)
        Hs = np.max(np.array(gn.f_Ham_terms(*q)) + 0*q[1], axis=0)
        Ms = np.max(np.array(gn.f_Mom_terms(*q)) + 0*q[1], axis=0)
        return np.max(np.abs(H)), np.max(np.abs(H))/np.max(Hs), np.max(np.abs(M)), np.max(np.abs(M))/np.max(Ms)

    # ---------- observer clock ----------
    def alh_at(self, y, xr, need_x=False):
        """alh (and d_x alh) at arbitrary points xr from the spectral state."""
        Ah, Bh, pA, pB = self.unpack(y)
        E = self.interp_setup(np.atleast_1d(xr))
        fA, fB = np.fft.rfft(Ah), np.fft.rfft(Bh)
        Ahr, Bhr = (E @ fA).real, (E @ fB).real
        z0, z1 = self.zeta_derivs(np.atleast_1d(xr), 1)
        alh = np.exp((np.log(Ahr) + 2*np.log(Bhr))/3 - z0)
        if not need_x:
            return alh
        Axr, Bxr = (E @ (fA*self.ik)).real, (E @ (fB*self.ik)).real
        return alh, alh*(-z1 + (Axr/Ahr + 2*Bxr/Bhr)/3)

    def rhs_forward(self, t, Y):
        """fields + tau (proper time of observer x_o): d tau/dt = alpha(t, x_o) = t alh."""
        y = Y[:-1]
        dy = self.rhs_fields(t, y)
        dtau = t*self.alh_at(y, self.x_o)[0]
        return np.concatenate([dy, [dtau]])

    def integrate(self, rhs, Y0, t0, t1, dense=False, events=None):
        return solve_ivp(rhs, (t0, t1), Y0, method='DOP853', rtol=self.rtol, atol=self.atol,
                         dense_output=dense, events=events)

    def evolve_forward(self):
        """Evolve from t_i to the observation time t_o (fixed observer proper time). Stores state at t_o."""
        # initial event: T_phys(t_init, x_o) = T_i = 1/t_i  <=>  t alh(t, x_o) = t_i (series solution, Newton)
        f = lambda tt: tt*self.series_fields(tt, np.array([self.x_o]))[4][0] - self.t_i
        t_init = brentq(f, self.t_i*(1 - 0.5), self.t_i*(1 + 0.5), xtol=1e-16, rtol=1e-15)
        # proper time from t_init to t_i (Gauss-Legendre on the series solution)
        gx, gw = leggauss(12)
        tt = 0.5*(self.t_i - t_init)*gx + 0.5*(self.t_i + t_init)
        tau_pre = 0.5*(self.t_i - t_init)*np.sum(gw*np.array([s*self.series_fields(s, np.array([self.x_o]))[4][0] for s in tt]))
        self.t_init, self.tau_pre = t_init, tau_pre
        self.tau_target = (self.t_obar**2 - self.t_i**2)/2 - tau_pre   # tau measured from t_i
        Y0 = np.concatenate([self.initial_state(), [0.0]])
        ev = lambda t, Y: Y[-1] - self.tau_target
        ev.terminal = True; ev.direction = 1
        sol = self.integrate(self.rhs_forward, Y0, self.t_i, self.t_obar*1.5, dense=False, events=ev)
        assert sol.status == 1, sol.message
        # polish t_o: re-integrate from last full step before the event, Newton on tau(t_o) - target
        ta, Ya = sol.t[-2], sol.y[:, -2]
        t_o = sol.t_events[0][0]
        for it in range(4):
            s2 = self.integrate(self.rhs_forward, Ya, ta, t_o)
            Yo = s2.y[:, -1]
            res = Yo[-1] - self.tau_target
            dres = t_o*self.alh_at(Yo[:-1], self.x_o)[0]
            t_o = t_o - res/dres
            if abs(res/dres) < 1e-15*t_o:
                break
        s2 = self.integrate(self.rhs_forward, Ya, ta, t_o)
        self.t_o, self.y_o = t_o, s2.y[:-1, -1]
        self.sol_forward = sol
        self.nsteps_forward = len(sol.t)
        self.T_o = 1.0/(t_o*self.alh_at(self.y_o, self.x_o)[0])
        if self.verbose:
            c = self.constraints(t_o, self.y_o)
            print(f"[eps={self.eps:+.4f}] forward: {len(sol.t)} steps, t_init={t_init:.12e}, t_o={t_o:.13f} "
                  f"(t_obar={self.t_obar}), T_o={self.T_o:.13e}; constraints at t_o: |H|={c[0]:.2e} rel={c[1]:.2e} |M|={c[2]:.2e} rel={c[3]:.2e}")
        return self

    # ---------- ray tracing ----------
    def ray_setup(self):
        """Observer tetrad: u = (1/alpha,0,0,0), e_x = d_x/A, e_y = d_y/B. Photon from sky direction n,
        p = E_o (u - n^a e_a), E_o = 1.  n = (mu, sqrt(1-mu^2), 0) in the triad (axisymmetry about x)."""
        mu, w = leggauss(self.nmu)
        self.mu, self.wmu = mu, w
        Ah, Bh, pA, pB = self.unpack(self.y_o)
        E = self.interp_setup(np.array([self.x_o]))
        Aho = (E @ np.fft.rfft(Ah)).real[0]; Bho = (E @ np.fft.rfft(Bh)).real[0]
        alho = self.alh_at(self.y_o, self.x_o)[0]
        A_o, B_o, al_o = self.t_o*Aho, self.t_o*Bho, self.t_o*alho
        px = -mu*A_o                       # p_x = A^2 p^x, p^x = -mu/A
        py = -np.sqrt(1 - mu**2)*B_o       # conserved
        pt = -al_o*np.ones(self.nmu)       # p_t = -alpha E
        self.py = py
        return np.concatenate([self.y_o, self.x_o + 0*mu, px, pt, 0*mu])

    def ray_quantities(self, t, Y):
        N, nr = self.N, self.nmu
        y = Y[:4*N]
        xr, px, pt = Y[4*N:4*N + nr], Y[4*N + nr:4*N + 2*nr], Y[4*N + 2*nr:4*N + 3*nr]
        Ah, Bh, pA, pB = self.unpack(y)
        E = self.interp_setup(xr)
        fA, fB, fpA, fpB = np.fft.rfft(Ah), np.fft.rfft(Bh), np.fft.rfft(pA), np.fft.rfft(pB)
        Ahr, Bhr = (E @ fA).real, (E @ fB).real
        Axr, Bxr = (E @ (fA*self.ik)).real, (E @ (fB*self.ik)).real
        pAr, pBr = (E @ fpA).real, (E @ fpB).real
        z0, z1 = self.zeta_derivs(xr, 1)
        alh = np.exp((np.log(Ahr) + 2*np.log(Bhr))/3 - z0)
        alh_x = alh*(-z1 + (Axr/Ahr + 2*Bxr/Bhr)/3)
        alh_t = alh**2*(pAr + 2*pBr)/3
        return y, xr, px, pt, Ahr, Bhr, Axr, Bxr, pAr, pBr, alh, alh_x, alh_t

    def rhs_rays(self, t, Y):
        """fields + rays, parametrised by coordinate time t (integrated backwards)."""
        y, xr, px, pt, Ahr, Bhr, Axr, Bxr, pAr, pBr, alh, alh_x, alh_t = self.ray_quantities(t, Y)
        dy = self.rhs_fields(t, y)
        py = self.py
        A, B, al = t*Ahr, t*Bhr, t*alh
        A_x, B_x, al_x = t*Axr, t*Bxr, t*alh_x
        A_t, B_t, al_t = Ahr + t*alh*Ahr*pAr, Bhr + t*alh*Bhr*pBr, alh + t*alh_t
        dt_dl = -pt/al**2
        dx_dl = px/A**2
        dpx_dl = -(pt**2*al_x/al**3 - px**2*A_x/A**3 - py**2*B_x/B**3)
        dpt_dl = -(pt**2*al_t/al**3 - px**2*A_t/A**3 - py**2*B_t/B**3)
        Eng = -pt/al
        c = (px/A)/Eng                                   # e_x component of propagation direction
        dS_dt = alh*(pAr - pBr)*(c**2 - 1/3)             # sigma_ee * dl/dt  (dl = E dlambda)
        return np.concatenate([dy, dx_dl/dt_dl, dpx_dl/dt_dl, dpt_dl/dt_dl, dS_dt])

    def emission_residual(self, t, Y):
        """alpha(t, x_r) - t_star  (= 0 where T_phys = T_*)."""
        q = self.ray_quantities(t, Y)
        return t*q[10] - self.t_star

    def emission_residual_dt(self, t, Y):
        y, xr, px, pt, Ahr, Bhr, Axr, Bxr, pAr, pBr, alh, alh_x, alh_t = self.ray_quantities(t, Y)
        dY = self.rhs_rays(t, Y)
        dxdt = dY[4*self.N:4*self.N + self.nmu]
        return alh + t*(alh_t + alh_x*dxdt)

    def trace_rays(self):
        Y0 = self.ray_setup()
        N, nr = self.N, self.nmu
        margin = 0.05 + 4*abs(self.eps)
        t_end = self.t_star*(1 - margin)
        sol = self.integrate(self.rhs_rays, Y0, self.t_o, t_end, dense=True)
        assert sol.status == 0, sol.message
        self.nsteps_back = len(sol.t)
        # locate emission events per ray
        res_steps = np.array([self.emission_residual(tt, sol.y[:, i]) for i, tt in enumerate(sol.t)])  # (nsteps, nr)
        t_e = np.zeros(nr); x_e = np.zeros(nr); px_e = np.zeros(nr); pt_e = np.zeros(nr); S_e = np.zeros(nr)
        E_e = np.zeros(nr); E_e_null = np.zeros(nr); alpha_e = np.zeros(nr); newton_res = np.zeros(nr)
        for r in range(nr):
            idx = np.where(np.sign(res_steps[:-1, r]) != np.sign(res_steps[1:, r]))[0]
            assert len(idx) >= 1, "no emission crossing for ray %d" % r
            i = idx[0]
            ta, tb = sol.t[i], sol.t[i + 1]
            g = lambda tt: self.emission_residual(tt, sol.sol(tt))[r]
            te = brentq(g, ta, tb, xtol=1e-14, rtol=1e-14)
            Ya = sol.y[:, i]
            for it in range(5):
                s2 = self.integrate(self.rhs_rays, Ya, ta, te)
                Ye = s2.y[:, -1]
                resid = self.emission_residual(te, Ye)[r]
                dres = self.emission_residual_dt(te, Ye)[r]
                step = -resid/dres
                te = te + step
                if abs(step) < 1e-15*abs(te):
                    break
            s2 = self.integrate(self.rhs_rays, Ya, ta, te)
            Ye = s2.y[:, -1]
            newton_res[r] = self.emission_residual(te, Ye)[r]
            q = self.ray_quantities(te, Ye)
            t_e[r] = te; x_e[r] = q[1][r]; px_e[r] = q[2][r]; pt_e[r] = q[3][r]
            S_e[r] = -Ye[4*N + 3*nr + r]        # integrated backwards in t: S = int_e^o sigma_ee dl = -(stored)
            al_e = te*q[10][r]; A_e = te*q[4][r]; B_e = te*q[5][r]
            alpha_e[r] = al_e
            E_e[r] = -pt_e[r]/al_e
            E_e_null[r] = np.sqrt(px_e[r]**2/A_e**2 + self.py[r]**2/B_e**2)
        self.t_e, self.x_e, self.S, self.E_e, self.E_e_null, self.newton_res = t_e, x_e, S_e, E_e, E_e_null, newton_res
        self.T_obs = (1.0/self.t_star)/E_e              # T_* E_o/E_e with E_o = 1
        self.T_obs_null = (1.0/self.t_star)/E_e_null
        self.ratio = self.T_obs/self.T_o                # should equal exp(-S) exactly (SPEC sec.3)
        self.monopole = 0.5*np.sum(self.wmu*self.T_obs)
        self.Theta = self.T_obs/self.monopole - 1
        self.exact_formula_residual = np.max(np.abs(self.ratio - np.exp(-self.S)))
        if self.verbose:
            print(f"[eps={self.eps:+.4f}] rays: {self.nsteps_back} steps back; |alpha_e - t_*|max={np.max(np.abs(alpha_e - self.t_star)):.1e}; "
                  f"null-cond E mismatch max={np.max(np.abs(E_e/E_e_null - 1)):.1e}; max|T_obs/T_o - exp(-S)|={self.exact_formula_residual:.2e}; "
                  f"Theta rms={np.sqrt(0.5*np.sum(self.wmu*self.Theta**2)):.3e}")
        return self

    def alm(self, f, lmax=30):
        """Legendre coefficients a_l0 = sqrt(4 pi (2l+1)) * (1/2) int f(mu) P_l(mu) dmu (real Y_l0 normalisation, axis = x)."""
        from scipy.special import eval_legendre
        return np.array([np.sqrt(4*np.pi*(2*l + 1))*0.5*np.sum(self.wmu*f*eval_legendre(l, self.mu)) for l in range(lmax + 1)])

    def run(self):
        return self.evolve_forward().trace_rays()
