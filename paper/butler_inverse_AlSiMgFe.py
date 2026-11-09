"""Inverse Butler: interface composition X_i^S(grad T) of the FCC_A1 nucleus of Al-0.8Si-0.6Mg-0.2Fe with the
interfacial energy imposed by the nucleation solution, gamma = gamma_SL(grad T), at T = T_L - Delta T(grad T):
    gamma_SL = gamma_i + RT/A_i ln(X_i^S/X_i^B) + (1/A_i)[G_i^{E,S}(T,X^S) - G_i^{E,B}(T,X^B)] - c X_i^S
Four equations plus sum X_i^S = 1 for the four X_i^S: the closure needs one more unknown, the configurational
free energy actually taken by the interface, c(grad T) (same X_i^S weighting as Eq. 47). c is found so that Butler
returns gamma_SL; its ratio to the nucleus value, c/Delta G_C, is the fraction of the substrate term spent at the
interface. Same thermodynamic inputs as butler_AlSiMgFe.py.
Writes outputs_AlSiMgFe/butler_inverse_AlSiMgFe.csv and Fig12_butler_surface_composition.png."""
import csv
import os

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import butler_AlSiMgFe as BU

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(BU.OUT, 'butler_inverse_AlSiMgFe.csv')
G_YIELD = 925.0          # end of the elastic range of the interface [K/m]
G_MEAN = 7902.38         # mean experimental thermal gradient [K/m]
NAMES = {'AL': 'Al', 'SI': 'Si', 'MG': 'Mg', 'FE': 'Fe'}
COL = {'AL': '#2a78d6', 'SI': '#eb6834', 'MG': '#1a9e5c', 'FE': '#8a3ffc'}
plt.rcParams.update({'font.size': 11, 'font.family': 'serif', 'mathtext.fontset': 'dejavuserif'})


def inverse(T, xb, gamma, gC, guess):
    """c in [Delta G_C, 0] such that Butler gives gamma; returns c, X^S and the solver state."""
    state = {'u': guess}

    def f(c):
        xs, g, u = BU.butler(T, xb, c, state['u'])
        state['u'], state['xs'] = u, xs
        return g - gamma

    c = brentq(f, gC, 0.0, xtol=1e-10, rtol=1e-10)
    xs, g, u = BU.butler(T, xb, c, state['u'])
    return c, xs, u


def solve():
    xb = BU.bulk_fractions()
    out = []
    for surf, pts in BU.load_sweep().items():
        u = None
        for p in pts:
            T = BU.T_L - p['DT']
            c, xs, u = inverse(T, xb, p['gamma'], p['G_C'], u)
            out.append([surf, p['grad'], p['DT'], T, p['gamma'], p['G_C'], c, c / p['G_C'], *xs])
        print(surf, 'done', len(pts))
    head = ['surface', 'grad', 'DT', 'T', 'gamma_nuc', 'G_C', 'c', 'c_over_GC'] + [f'XS_{e}' for e in BU.EL]
    with open(CSV, 'w', newline='') as f:
        w = csv.writer(f, delimiter=';')
        w.writerow(head)
        w.writerows([[r[0]] + [f'{v:.6e}' for v in r[1:]] for r in out])
    print('written', CSV)


def load():
    data = {}
    with open(CSV, newline='') as f:
        for r in csv.DictReader(f, delimiter=';'):
            d = data.setdefault(r['surface'], {})
            for k, v in r.items():
                if k != 'surface':
                    d.setdefault(k, []).append(float(v))
    return {s: {k: np.asarray(v) for k, v in d.items()} for s, d in data.items()}


def figure(data):
    e, p = data['elastic'], data['elasto-plastic']
    xs_std = BU.butler(BU.T_L, BU.bulk_fractions())[0]
    fig, axs = plt.subplots(1, 2, figsize=(11.0, 4.5))

    ax = axs[0]
    for k, x in enumerate(BU.EL):
        ax.plot(e['grad'], e[f'XS_{x}'], color=COL[x], lw=2, label=f'{NAMES[x]}')
        ax.plot(p['grad'], p[f'XS_{x}'], color=COL[x], lw=1.2, ls='--')
        ax.axhline(xs_std[k], color=COL[x], lw=1, ls=':')
    ax.plot([], [], color='0.3', lw=2, label='elastic')
    ax.plot([], [], color='0.3', lw=1.2, ls='--', label='elasto-plastic')
    ax.plot([], [], color='0.3', lw=1, ls=':', label='standard Butler')
    ax.set_ylabel(r'Interface composition $X_i^S$ [-]')
    ax.set_yscale('log')
    ax.set_ylim(1e-5, 1.5)

    ax = axs[1]
    ax.plot(e['grad'], e['c_over_GC'], color='#2a78d6', lw=2, label='elastic')
    ax.plot(p['grad'], p['c_over_GC'], color='#eb6834', lw=2, ls='--', label='elasto-plastic')
    ax.set_ylabel(r'$c/\Delta G_C$ [-]')
    ax.set_ylim(0, 1)

    for k, ax in enumerate(axs):
        ax.axvspan(e['grad'][0], G_YIELD, color='#2a78d6', alpha=0.08, lw=0)
        ax.axvline(G_YIELD, color='k', lw=1, ls='-.', label='Yield: elastic range (shaded) ends' if k == 1 else None)
        ax.axvline(G_MEAN, color='0.4', lw=1, ls=':', label=r'Mean experimental $\nabla T$' if k == 1 else None)
        ax.set_xscale('log')
        ax.set_xlabel(r'Thermal gradient $\nabla T$ [K m$^{-1}$]')
        ax.grid(True, ls='--', color='0.8', lw=0.5)
        ax.text(0.02, 0.97, f'({"ab"[k]})', transform=ax.transAxes, va='top', fontweight='bold')
        ax.legend(fontsize=7.5, loc='lower right', ncol=2 if k == 0 else 1, frameon=True, edgecolor='k')
    fig.tight_layout()
    path = os.path.join(HERE, 'Fig12_butler_surface_composition.png')
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print('saved', path)
    for surf, d in data.items():
        for g in (d['grad'][0], G_YIELD, G_MEAN, 1e5, d['grad'][-1]):
            vals = {k: float(np.interp(np.log(g), np.log(d['grad']), d[k]))
                    for k in ['DT', 'gamma_nuc', 'c_over_GC'] + [f'XS_{x}' for x in BU.EL]}
            print(f'{surf:15s} gradT={g:10.4g} ' + ' '.join(f'{k}={v:.4g}' for k, v in vals.items()))


if __name__ == '__main__':
    import sys
    if '--plot' not in sys.argv:
        solve()
    figure(load())
