"""sphere.py -- Gauss-Legendre x uniform-phi quadrature on the sky and spherical-harmonic projection (Condon-Shortley, scipy)."""
import numpy as np
from numpy.polynomial import legendre as Lg
from scipy.special import sph_harm_y

class Sky:
    def __init__(self, nth=48, nph=None):
        nph = nph or 2 * nth
        x, w = Lg.leggauss(nth)
        th = np.arccos(x); ph = 2 * np.pi * np.arange(nph) / nph
        TH, PH = np.meshgrid(th, ph, indexing='ij')
        self.th, self.ph = TH.ravel(), PH.ravel()
        self.w = (np.repeat(w, nph) * 2 * np.pi / nph)
        self.nv = np.stack([np.sin(self.th) * np.cos(self.ph), np.sin(self.th) * np.sin(self.ph), np.cos(self.th)], axis=1)
        self._Y = {}
    def Y(self, l, m):
        if (l, m) not in self._Y:
            self._Y[(l, m)] = sph_harm_y(l, m, self.th, self.ph)
        return self._Y[(l, m)]
    def alm(self, f, lmax):
        """a_lm = int dOmega Y*_lm f  -> dict l -> array over m=-l..l."""
        return {l: np.array([np.sum(self.w * np.conj(self.Y(l, m)) * f) for m in range(-l, l + 1)]) for l in range(lmax + 1)}
    def power(self, f, lmax):
        a = self.alm(f, lmax)
        return np.array([np.sum(np.abs(a[l])**2) for l in range(lmax + 1)])
