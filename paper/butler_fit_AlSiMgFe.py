"""Fit of Butler's equation to the nucleation solution: which parameter makes gamma_Butler = gamma_SL(grad T)?

General configurational term (both series of Section 5.9):
    G_i^{E,C}/A_i = k [ (1 - w)/n + w X_i^S ] Delta G_C^x,   x = conf (Delta G_C) or rel (Delta G_C (1 - gamma_0/gamma))
    w = 1: weighting by the interface fraction (Figures 11-12);  w = 0: the same term for every component.
Fits (least squares of gamma_Butler - gamma_SL over the elastic sweep):
    A  the pure-component energies gamma_i (4 parameters), w = 1, k = 1
    B  w, k and a uniform shift delta of the gamma_i (3 parameters)
Each fit is checked on the whole elastic and elasto-plastic sweeps.
"""
import numpy as np
from scipy.optimize import least_squares

import butler_AlSiMgFe as BU

XB = BU.bulk_fractions()
SW = BU.load_sweep()
N_FIT = 10


def solve(T, gC, gam_i, w, k):
    """Butler with gamma_i replaced and the general configurational term; returns gamma and X^S."""
    saved = dict(BU.GAMMA_I)
    BU.GAMMA_I.update(dict(zip(BU.EL, gam_i)))
    try:
        if w == 1.0:
            xs, g, _ = BU.butler(T, XB, k * gC)
        else:
            # uniform part shifts all equations equally: gamma -> gamma - k (1-w)/n gC, X^S unchanged by it
            xs, g, _ = BU.butler(T, XB, w * k * gC) if w > 0 else BU.butler(T, XB, 0.0)
            g -= k * (1 - w) / len(BU.EL) * gC
    finally:
        BU.GAMMA_I.clear()
        BU.GAMMA_I.update(saved)
    return g, xs


def gc_of(p, series):
    return p['G_C'] if series == 'conf' else p['G_C'] * (1 - BU.GAMMA_0 / p['gamma'])


def residuals(pts, series, gam_i, w, k):
    out = []
    for p in pts:
        try:
            out.append(solve(BU.T_L - p['DT'], gc_of(p, series), gam_i, w, k)[0] - p['gamma'])
        except (RuntimeError, ValueError):
            out.append(10.0)
    return np.array(out)


def report(tag, series, gam_i, w, k):
    print(f'  {tag}: gamma_i = {np.round(gam_i, 4)}, w = {w:.4f}, k = {k:.4f}')
    for surf in ('elastic', 'elasto-plastic'):
        r = residuals(SW[surf], series, gam_i, w, k)
        print(f'    {surf:15s} max |gamma_Butler - gamma_SL| = {np.max(np.abs(r)):.3e} J/m2, '
              f'at start {r[0]:+.4f}, at end {r[-1]:+.4f}')


def main():
    pts = SW['elastic']
    idx = np.unique(np.geomspace(1, len(pts), N_FIT).astype(int) - 1)
    fit_pts = [pts[i] for i in idx]
    g0 = np.array([BU.GAMMA_I[e] for e in BU.EL])
    print('Turnbull gamma_i:', np.round(g0, 4))
    for series in ('conf', 'rel'):
        print(f'\n=== series {series} ===')
        report('as in Figures 11-12', series, g0, 1.0, 1.0)

        fa = least_squares(lambda q: residuals(fit_pts, series, q, 1.0, 1.0), g0, bounds=(-5, 20), x_scale=0.1)
        report('fit A (gamma_i, w = 1)', series, fa.x, 1.0, 1.0)

        def rb(q):
            return residuals(fit_pts, series, g0 + q[2], q[0], q[1])
        fb = least_squares(rb, [0.5, 0.5, 0.0], bounds=([0, 0, -1], [1, 2, 1]), x_scale=[0.1, 0.1, 0.01])
        report(f'fit B (w, k, delta = {fb.x[2]:+.4f})', series, g0 + fb.x[2], fb.x[0], fb.x[1])


if __name__ == '__main__':
    main()
