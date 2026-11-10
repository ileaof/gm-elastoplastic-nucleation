"""Butler's formulation for the thermodynamic equilibrium shift (Ferreira, Materials 2026, 19, 4676, Eq. 34) with the
excess configurational free energy computed as a genuine excess function of the surface solution: the atomic-size
mismatch (hard-sphere) contribution of Mansoori, Carnahan, Starling and Leland (J. Chem. Phys. 54, 1523, 1971) - the
size-mismatch mechanism that Section 3.2 of the article invokes for the configurational free energy (Vitreloy 1).

    gamma = gamma_i + RT/A_i ln(X_i^S/X_i^B) + (1/A_i)[G_i^{E,S}(T,X^S) - G_i^{E,B}(T,X^B) - G_i^{E,C}(T,X^S)],  sum X^S = 1
    G^{E,C}(T,x) = RT [A^ex_MCSL(x, xi_S) - A^ex_pure(xi_S)]   (J/mol; zero for equal atomic sizes)
    G_i^{E,C}    = partial molar quantity of G^{E,C} (same routine as G^{E,S}, G^{E,B})

Forward Butler: gamma and X^S are the unknowns (Butler and nucleation are different things); the thermal gradient enters
through the temperature of nucleation, T = T_L - Delta T(grad T), of the elastic and elasto-plastic solutions.
Surface packing fraction xi_S = 0.60 (between the liquid, ~0.45, and FCC, 0.74); sensitivity 0.45-0.74 at T_L.
Atomic diameters: twice the metallic radii for coordination 12.
Liquid excess (COST507), beta, Turnbull gamma_i and A_i as in butler_composition_ligas.py.

Writes outputs_AlSiMgFe/butler_configurational_<alloy>.csv and Fig_butler_configurational_<alloy>.png.
"""
import csv
import os
import sys
from math import log

import numpy as np
from scipy.optimize import brentq, fsolve

from butler_composition_ligas import ALLOYS, BETA, OUT, R, Alloy, read_liquid_parameters

HERE = os.path.dirname(os.path.abspath(__file__))
R_METAL = {'AL': 1.432, 'SI': 1.319, 'MG': 1.602, 'FE': 1.274, 'CU': 1.278, 'NB': 1.468}   # Angstrom, CN 12
XI_S = 0.60
XI_SENS = (0.45, 0.52, 0.60, 0.67, 0.74)


def a_ex_mcsl(x, d, xi):
    """Excess Helmholtz energy of the hard-sphere mixture over the ideal gas, per particle / kT (MCSL 1971)."""
    x = np.clip(np.asarray(x, float), 0.0, None)
    x = x / x.sum()
    s3 = float(x @ d ** 3)
    xr = x * d ** 3 / s3                                     # volume fractions xi_i/xi
    y1 = y2 = 0.0
    n = len(x)
    for i in range(n):
        for j in range(i + 1, n):
            dij = np.sqrt(xr[i] * xr[j]) * (d[i] - d[j]) ** 2 / (d[i] * d[j]) * np.sqrt(x[i] * x[j])
            y1 += dij * (d[i] + d[j]) / np.sqrt(d[i] * d[j])
            y2 += dij * sum(xr[k] * np.sqrt(d[i] * d[j]) / d[k] for k in range(n))
    y3 = float(sum(xr[i] ** (2 / 3) * x[i] ** (1 / 3) for i in range(n))) ** 3
    z = 1.0 - xi
    return (-1.5 * (1 - y1 + y2 + y3) + (3 * y2 + 2 * y3) / z + 1.5 * (1 - y1 - y2 - y3 / 3) / z ** 2
            + (y3 - 1) * log(z))


def a_ex_pure(xi):
    z = 1.0 - xi
    return -3.0 + 2.0 / z + 1.0 / z ** 2                   # Carnahan-Starling: (4 xi - 3 xi^2)/(1 - xi)^2


class Configurational:
    def __init__(self, els, xi):
        self.d = np.array([2 * R_METAL[e] for e in els])
        self.xi = xi

    def g(self, x, T):
        return R * T * (a_ex_mcsl(x, self.d, self.xi) - a_ex_pure(self.xi))

    def partial(self, x, T, h=1e-7):
        grad = np.empty(len(x))
        for k in range(len(x)):
            xp, xm = x.copy(), x.copy()
            xp[k] += h
            xm[k] -= h
            grad[k] = (self.g(xp, T) - self.g(xm, T)) / (2 * h)
        return self.g(x, T) + grad - x @ grad


def butler_forward(al, T, conf=None, guess=None):
    """Butler with (conf) or without the configurational term: gamma and X^S at temperature T."""
    xb, A, g0, liq = al.xb, al.A, al.g0, al.liq
    gb = liq.partial(xb, T)
    xs = xb.copy() if guess is None else guess.copy()
    for _ in range(300 if guess is None else 60):           # damped fixed point on the surface terms
        gs = BETA * liq.partial(xs, T) - (conf.partial(xs, T) if conf else 0.0)
        xk = lambda g: xb * np.exp(A * (g - g0 - (gs - gb) / A) / (R * T))
        gam = brentq(lambda g: xk(g).sum() - 1.0, -20.0, 20.0, xtol=1e-15)
        xs = 0.5 * xs + 0.5 * xk(gam) / xk(gam).sum()

    def res(u):
        x = np.concatenate(([1 - np.exp(u[:-1]).sum()], np.exp(u[:-1])))
        if x[0] <= 0:
            return np.full(len(u), 1e3)
        gsx = BETA * liq.partial(x, T) - (conf.partial(x, T) if conf else 0.0)
        return g0 + R * T / A * np.log(x / xb) + (gsx - gb) / A - u[-1]

    u, _, ier, msg = fsolve(res, np.concatenate((np.log(xs[1:]), [gam])), full_output=True, xtol=1e-14)
    if np.max(np.abs(res(u))) > 1e-8:
        raise RuntimeError(f'{al.key}: Butler did not converge at T={T:.3f} K: {msg}')
    x = np.concatenate(([1 - np.exp(u[:-1]).sum()], np.exp(u[:-1])))
    gc_area = (conf.partial(x, T) / A) if conf else np.zeros(len(x))
    return u[-1], x, gc_area


def run(key, pars):
    al = Alloy(key, pars)
    conf = Configurational(al.els, XI_S)
    g_std, x_std, _ = butler_forward(al, al.T_L)
    g_c0, x_c0, gc0 = butler_forward(al, al.T_L, conf)
    print(f'\n== {al.title}  T_L = {al.T_L} K   xi_S = {XI_S}  d [A] =', dict(zip(al.els, conf.d)))
    print(f'   at T_L: standard Butler gamma = {g_std:.5f} J/m2; with G^E,C gamma = {g_c0:.5f} J/m2')
    print('      X^S standard :', dict(zip(al.els, np.round(x_std, 6))))
    print('      X^S with G^EC:', dict(zip(al.els, np.round(x_c0, 6))))
    print('      G_i^E,C/A_i [J/m2]:', dict(zip(al.els, np.round(gc0, 5))),
          f' G^E,C(X^S) = {conf.g(x_c0, al.T_L):.1f} J/mol')
    sens = []
    for xi in XI_SENS:
        g, x, _ = butler_forward(al, al.T_L, Configurational(al.els, xi), x_c0)
        sens.append((xi, g, x))
        print(f'      xi_S = {xi:.2f}: gamma = {g:.5f} J/m2, X^S =', dict(zip(al.els, np.round(x, 5))))

    out = []
    for surf, pts in al.load().items():
        gs_, gc_ = x_std, x_c0
        for grad, DT, gam_nuc, _ in pts:
            T = al.T_L - DT
            g0_, xs0, _ = butler_forward(al, T, None, gs_)
            g1_, xs1, gca = butler_forward(al, T, conf, gc_)
            gs_, gc_ = xs0, xs1
            out.append([surf, grad, DT, T, gam_nuc, g0_, g1_, *xs0, *xs1, *gca])
    path = os.path.join(OUT, f'butler_configurational_{key}.csv')
    with open(path, 'w', newline='') as f:
        w = csv.writer(f, delimiter=';')
        w.writerow(['surface', 'grad', 'DT', 'T', 'gamma_nuc', 'gamma_std', 'gamma_conf'] +
                   [f'XS_std_{e}' for e in al.els] + [f'XS_conf_{e}' for e in al.els] +
                   [f'GEC_over_A_{e}' for e in al.els])
        w.writerows([[r[0]] + [f'{v:.6e}' for v in r[1:]] for r in out])
    n = len(al.els)
    for surf in dict.fromkeys(r[0] for r in out):
        sub = [r for r in out if r[0] == surf]
        for r in (sub[0], sub[len(sub) // 2], sub[-1]):
            print(f'   {surf:15s} gradT={r[1]:9.3g} T={r[3]:8.3f} gamma std/conf = {r[5]:.5f}/{r[6]:.5f}  X^S conf ' +
                  ' '.join(f'{e}={v:.3e}' for e, v in zip(al.els, r[7 + n:7 + 2 * n])))
    figure(al, out, sens)


def figure(al, out, sens):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'serif', 'font.size': 11, 'mathtext.fontset': 'dejavuserif'})
    col = {'AL': '#2a7ae2', 'SI': '#eb6834', 'MG': '#1a9e5c', 'FE': '#8a3ffc', 'CU': '#c0392b', 'NB': '#7f6000'}
    nm = {'AL': 'Al (SV)', 'SI': 'Si', 'MG': 'Mg', 'FE': 'Fe', 'CU': 'Cu', 'NB': 'Nb'}
    n = len(al.els)
    surfs = list(dict.fromkeys(r[0] for r in out))
    style = {'elastic': '-', 'elasto-plastic': '--'}
    fig, ax = plt.subplots(1, 3, figsize=(17.5, 5.0))
    for surf in surfs:
        sub = [r for r in out if r[0] == surf]
        g = [r[1] for r in sub]
        for k, e in enumerate(al.els):
            ax[0].plot(g, [r[7 + n + k] for r in sub], style[surf], color=col[e], lw=2.0 if surf == surfs[0] else 1.4,
                       label=nm[e] if surf == surfs[0] else None)
            ax[0].plot(g, [r[7 + k] for r in sub], ':', color=col[e], lw=1.0)
            ax[2].plot(g, [r[7 + 2 * n + k] for r in sub], style[surf], color=col[e], lw=1.8,
                       label=nm[e] if surf == surfs[0] else None)
        ax[1].plot(g, [r[6] for r in sub], style[surf], color='#c0392b', lw=2.0, label=f'with $G^{{E,C}}$, {surf}')
        ax[1].plot(g, [r[5] for r in sub], ':', color='0.3', lw=1.4, label=f'standard Butler, {surf}')
    for surf in surfs:
        ax[0].plot([], [], 'k' + style[surf], label=f'with $G^{{E,C}}$, {surf}')
    ax[0].plot([], [], 'k:', label='standard Butler')
    ax[0].set_yscale('log')
    ax[0].set_ylabel('Interface composition $X_i^S$ [-]')
    ax[1].set_ylabel(r'Interfacial energy (Butler) $\gamma$ [J m$^{-2}$]')
    ax[2].set_ylabel(r'Configurational term $G_i^{E,C}/A_i$ [J m$^{-2}$]')
    ax[2].axhline(0, color='0.5', lw=0.6)
    txt = '\n'.join(f'$\\xi_S$={xi:.2f}: $\\gamma$={g:.4f}' for xi, g, _ in sens)
    ax[1].text(0.97, 0.05, 'at $T_L$:\n' + txt, transform=ax[1].transAxes, ha='right', va='bottom', fontsize=8,
               bbox=dict(facecolor='white', edgecolor='0.7'))
    for a, t in zip(ax, ('(a)', '(b)', '(c)')):
        a.set_xscale('log')
        a.set_xlabel(r'Thermal gradient $\nabla T$ [K m$^{-1}$]')
        a.grid(True, which='both', ls='--', lw=0.4, alpha=0.6)
        a.text(0.03, 0.95, t, transform=a.transAxes, fontweight='bold', va='top')
        a.legend(fontsize=8, loc='upper right' if a is not ax[1] else 'upper left')
    fig.suptitle(f'{al.title} (FCC_A1): Butler with the equilibrium shift, excess configurational free energy from the '
                 f'atomic-size mismatch (hard spheres, $\\xi_S$ = {XI_S}); $T = T_L - \\Delta T(\\nabla T)$', fontsize=11)
    fig.tight_layout()
    path = os.path.join(HERE, f'Fig_butler_configurational_{al.key}.png')
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print('   written', path)


if __name__ == '__main__':
    pars = read_liquid_parameters()
    for key in (sys.argv[1:] or list(ALLOYS)):
        run(key, pars)
