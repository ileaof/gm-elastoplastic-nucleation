"""Butler with an element-specific mechanical term: interface composition X_i^S(grad T) of the FCC_A1 nucleus of
Al-0.8Si-0.6Mg-0.2Fe driven by the size misfit of each element in the stressed Gurtin-Murdoch surface.
    gamma = gamma_i + delta + RT/A_i ln(X_i^S/X_i^B) + (1/A_i)[G_i^{E,S} - G_i^{E,B}] - Delta G_C^rel/n
            + phi_i - sum_j X_j^S phi_j,        phi_i = -S(grad T) eps_i,   n = 4,   T = T_L - Delta T
  - uniform configurational term Delta G_C^rel/n = Delta G_C (1 - gamma_0/gamma)/n (the same for every element) and
    delta = gamma_0 - gamma_std(T_L): Butler equals the nucleation gamma_SL (fit B of butler_fit_AlSiMgFe.py);
  - zero-mean mechanical term: the mechanical energy already in gamma_SL is not counted twice, only its partition
    between the elements changes X^S (surface analogue of the McLean/Eshelby interaction p dV);
  - areal misfit eps_i = (r_i/r_Al)^2 - 1 from the metallic radii for coordination 12;
  - S, the surface stress felt by the solute, three variants: 'tau' total surface tension, 'tau_rel' = tau - sigma_0
    (zero at equilibrium), 'elastic' = excess elastic stress s = tau_e + |sigma_0| (bounded by sigma_y when yielding).
Writes outputs_AlSiMgFe/butler_mech_AlSiMgFe.csv and Fig12_butler_mechanical_composition.png."""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import butler_AlSiMgFe as BU

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(BU.OUT, 'butler_mech_AlSiMgFe.csv')
SIGMA_0 = 1.09                                                   # |sigma_0| of the nucleation model [N/m]
R_CN12 = {'AL': 1.432, 'SI': 1.319, 'MG': 1.602, 'FE': 1.274}   # metallic radii, coordination 12 [Angstrom]
EPS = np.array([(R_CN12[e] / R_CN12['AL']) ** 2 - 1 for e in BU.EL])
VARIANTS = {'tau': r'total surface tension $\tau$', 'tau_rel': r'$\tau-\sigma_0$',
            'elastic': r'elastic excess stress $s$'}
G_YIELD, G_MEAN = 925.0, 7902.38
NAMES = {'AL': 'Al', 'SI': 'Si', 'MG': 'Mg', 'FE': 'Fe'}
COL = {'AL': '#2a78d6', 'SI': '#eb6834', 'MG': '#1a9e5c', 'FE': '#8a3ffc'}
plt.rcParams.update({'font.size': 11, 'font.family': 'serif', 'mathtext.fontset': 'dejavuserif'})


def load_sweep():
    rows = {}
    with open(BU.SWEEP, newline='') as f:
        for r in csv.DictReader(f, delimiter=';'):
            if r['G_C'] in ('', 'nan'):
                continue
            rows.setdefault(r['surface'], []).append(
                {k: float(r[k]) for k in ('grad', 'DT', 'gamma', 'G_C', 'tau', 'tau_e')})
    return rows


def stress(p, v):
    return {'tau': p['tau'], 'tau_rel': p['tau'] + SIGMA_0, 'elastic': p['tau_e'] + SIGMA_0}[v]


def solve():
    xb = BU.bulk_fractions()
    delta = BU.GAMMA_0 - BU.butler(BU.T_L, xb)[1]
    print('eps_i =', dict(zip(BU.EL, np.round(EPS, 4))), ' delta =', round(delta, 5))
    for e in BU.EL:
        BU.GAMMA_I[e] += delta
    out = []
    for surf, pts in load_sweep().items():
        for v in VARIANTS:
            u = None
            for p in pts:
                T = BU.T_L - p['DT']
                S = stress(p, v)
                xs, g, u = BU.butler(T, xb, 0.0, u, phi=-S * EPS)
                g -= p['G_C'] * (1 - BU.GAMMA_0 / p['gamma']) / len(BU.EL)
                out.append([surf, v, p['grad'], p['DT'], S, p['gamma'], g, *xs])
            print(surf, v, 'done')
    head = ['surface', 'variant', 'grad', 'DT', 'S', 'gamma_nuc', 'gamma_butler'] + [f'XS_{e}' for e in BU.EL]
    with open(CSV, 'w', newline='') as f:
        w = csv.writer(f, delimiter=';')
        w.writerow(head)
        w.writerows([r[:2] + [f'{x:.6e}' for x in r[2:]] for r in out])
    print('written', CSV)


def load():
    data = {}
    with open(CSV, newline='') as f:
        for r in csv.DictReader(f, delimiter=';'):
            d = data.setdefault((r['surface'], r['variant']), {})
            for k, x in r.items():
                if k not in ('surface', 'variant'):
                    d.setdefault(k, []).append(float(x))
    return {key: {k: np.asarray(x) for k, x in d.items()} for key, d in data.items()}


def figure(data):
    xs_std = BU.butler(BU.T_L, BU.bulk_fractions())[0]
    fig, axs = plt.subplots(1, 3, figsize=(15.0, 4.6), sharey=True)
    for k, (v, title) in enumerate(VARIANTS.items()):
        ax = axs[k]
        e, p = data[('elastic', v)], data[('elasto-plastic', v)]
        for j, x in enumerate(BU.EL):
            ax.plot(e['grad'], e[f'XS_{x}'], color=COL[x], lw=2, label=NAMES[x])
            ax.plot(p['grad'], p[f'XS_{x}'], color=COL[x], lw=1.2, ls='--')
            ax.axhline(xs_std[j], color=COL[x], lw=1, ls=':')
        ax.plot([], [], color='0.3', lw=2, label='elastic')
        ax.plot([], [], color='0.3', lw=1.2, ls='--', label='elasto-plastic')
        ax.plot([], [], color='0.3', lw=1, ls=':', label='standard Butler')
        ax.axvspan(e['grad'][0], G_YIELD, color='#2a78d6', alpha=0.08, lw=0)
        ax.axvline(G_YIELD, color='k', lw=1, ls='-.', label='Yield: elastic range (shaded) ends')
        ax.axvline(G_MEAN, color='0.4', lw=1, ls=':', label=r'Mean experimental $\nabla T$')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(1e-6, 1.5)
        ax.set_xlabel(r'Thermal gradient $\nabla T$ [K m$^{-1}$]')
        ax.set_title(f'S = {title}', fontsize=10.5)
        ax.grid(True, ls='--', color='0.8', lw=0.5)
        ax.text(0.02, 0.97, f'({"abc"[k]})', transform=ax.transAxes, va='top', fontweight='bold')
        if k == 0:
            ax.set_ylabel(r'Interface composition $X_i^S$ [-]')
            ax.legend(fontsize=7.5, loc='lower right', ncol=2, frameon=True, edgecolor='k')
    fig.tight_layout()
    path = os.path.join(HERE, 'Fig12_butler_mechanical_composition.png')
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print('saved', path)
    for (surf, v), d in sorted(data.items(), key=lambda q: (q[0][1], q[0][0])):
        print(f'{v:8s} {surf:15s} max |gamma_Butler - gamma_SL| = {np.max(np.abs(d["gamma_butler"] - d["gamma_nuc"])):.2e}')
        for g in (d['grad'][0], G_YIELD, G_MEAN, 1e5, d['grad'][-1]):
            vals = {k: float(np.interp(np.log(g), np.log(d['grad']), d[k])) for k in ['S'] + [f'XS_{x}' for x in BU.EL]}
            print(f'    gradT={g:10.4g} ' + ' '.join(f'{k}={x:.4g}' for k, x in vals.items()))


if __name__ == '__main__':
    import sys
    if '--plot' not in sys.argv:
        solve()
    figure(load())
