"""Figures 11-12 and Tables 11-12 of the paper: Butler equation of the FCC_A1 nucleus of Al-0.8Si-0.6Mg-0.2Fe solved along
the thermal gradient (standard Butler vs Butler with equilibrium shift, elastic and elasto-plastic surfaces).
Series 'shift': configurational term Delta G_C of the nucleus (Figure 11); series 'rel': the same term counted from
equilibrium, Delta G_C (1 - gamma_0/gamma), so that Butler starts at the start of nucleation (Figure 12).
Reads outputs_AlSiMgFe/butler_AlSiMgFe.csv written by butler_AlSiMgFe.py."""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, 'outputs_AlSiMgFe', 'butler_AlSiMgFe.csv')
plt.rcParams.update({'font.size': 11, 'font.family': 'serif', 'mathtext.fontset': 'dejavuserif'})

G_YIELD = 925.0          # end of the elastic range of the interface [K/m]
G_MEAN = 7902.38         # mean experimental thermal gradient [K/m]
EL = ('AL', 'SI', 'MG', 'FE')
NAMES = {'AL': 'Al', 'SI': 'Si', 'MG': 'Mg', 'FE': 'Fe'}
COL = {'AL': '#2a78d6', 'SI': '#eb6834', 'MG': '#1a9e5c', 'FE': '#8a3ffc'}


def load():
    data = {}
    with open(CSV, newline='') as f:
        for r in csv.DictReader(f, delimiter=';'):
            d = data.setdefault(r['surface'], {})
            for k, v in r.items():
                if k != 'surface':
                    d.setdefault(k, []).append(float(v))
    return {s: {k: np.asarray(v) for k, v in d.items()} for s, d in data.items()}


def at(d, key, grad):
    """Value of key at the thermal gradient grad (interpolation in log grad)."""
    return float(np.interp(np.log(grad), np.log(d['grad']), d[key]))


def table_rows(data, series='shift'):
    """Rows of Table 10: elastic / elasto-plastic at the yield point, the mean experimental gradient and 1e5 K/m."""
    e, p = data['elastic'], data['elasto-plastic']
    rows = []
    for g, name in ((G_YIELD, 'yield'), (G_MEAN, 'mean'), (1e5, '1e5')):
        def two(key, fmt):
            return f'{at(e, key, g):{fmt}} / {at(p, key, g):{fmt}}'
        rows.append((name, g, {
            'DT': two('DT', '.4f'), 'gamma_nuc': two('gamma_nuc', '.4f'),
            'std': two('gamma_butler_std', '.4f'), 'shift': two(f'gamma_butler_{series}', '.4f'),
            **{f'XS_{x}': two(f'XS_{series}_{x}', '.4f') for x in EL},
            'offset': two_diff(e, p, g, series)}))
    return rows


def two_diff(e, p, g, series='shift'):
    de = at(e, f'gamma_butler_{series}', g) - at(e, 'gamma_nuc', g)
    dp = at(p, f'gamma_butler_{series}', g) - at(p, 'gamma_nuc', g)
    return f'{de:.4f} / {dp:.4f}'


def figure(data, series, fname):
    e, p = data['elastic'], data['elasto-plastic']
    tag = 'with shift' if series == 'shift' else 'from equilibrium'
    fig, axs = plt.subplots(1, 3, figsize=(15.0, 4.5))

    ax = axs[0]
    ax.plot(e['grad'], e['gamma_nuc'], color='0.25', lw=1.6, label=r'Nucleation, elastic $\gamma_{SL}^{el}$')
    ax.plot(p['grad'], p['gamma_nuc'], color='0.25', lw=1.6, ls='--', label=r'Nucleation, elasto-plastic $\gamma_{SL}^{ep}$')
    ax.plot(e['grad'], e[f'gamma_butler_{series}'], color='#2a78d6', lw=2, label=f'Butler {tag}, elastic')
    ax.plot(p['grad'], p[f'gamma_butler_{series}'], color='#eb6834', lw=2, ls='--', label=f'Butler {tag}, elasto-plastic')
    ax.plot(e['grad'], e['gamma_butler_std'], color='#1a9e5c', lw=2, ls=':', label='Standard Butler (no shift)')
    ax.set_ylabel(r'Interfacial energy $\gamma$ [J m$^{-2}$]')
    ax.set_yscale('log')

    ax = axs[1]
    for x in EL:
        ax.plot(p['grad'], p[f'XS_{series}_{x}'], color=COL[x], lw=2, label=f'{NAMES[x]}, {tag} (elasto-plastic)')
        ax.plot(e['grad'], e[f'XS_{series}_{x}'], color=COL[x], lw=1, ls='--')
        ax.axhline(e[f'XS_std_{x}'][0], color=COL[x], lw=1, ls=':')
    ax.plot([], [], color='0.3', lw=1, ls='--', label=f'{tag} (elastic)')
    ax.plot([], [], color='0.3', lw=1, ls=':', label='standard Butler')
    ax.set_ylabel(r'Interface composition $X_i^S$ [-]')
    ax.set_yscale('log')
    ax.set_ylim(1e-5, 1.5)

    ax = axs[2]
    ax.plot(e['grad'], e[f'gamma_butler_{series}'] - e['gamma_nuc'], color='#2a78d6', lw=2, label='elastic')
    ax.plot(p['grad'], p[f'gamma_butler_{series}'] - p['gamma_nuc'], color='#eb6834', lw=2, ls='--', label='elasto-plastic')
    ax.axhline(e['gamma_butler_std'][0], color='#1a9e5c', lw=1.5, ls=':', label='standard Butler')
    ax.set_ylabel(r'$\gamma_{Butler}-\gamma_{SL}$ [J m$^{-2}$]')

    for k, ax in enumerate(axs):
        ax.axvspan(e['grad'][0], G_YIELD, color='#2a78d6', alpha=0.08, lw=0)
        ax.axvline(G_YIELD, color='k', lw=1, ls='-.', label='Yield: elastic range (shaded) ends' if k == 2 else None)
        ax.axvline(G_MEAN, color='0.4', lw=1, ls=':', label=r'Mean experimental $\nabla T$' if k == 2 else None)
        ax.set_xscale('log')
        ax.set_xlabel(r'Thermal gradient $\nabla T$ [K m$^{-1}$]')
        ax.grid(True, ls='--', color='0.8', lw=0.5)
        ax.text(0.02, 0.97, f'({"abc"[k]})', transform=ax.transAxes, va='top', fontweight='bold')
        ax.legend(fontsize=7.5, loc='lower right' if k != 1 else 'lower left', frameon=True, edgecolor='k')
    fig.tight_layout()
    path = os.path.join(HERE, fname)
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print('saved', path)
    for name, g, r in table_rows(data, series):
        print(series, name, g, r)


def main():
    data = load()
    figure(data, 'shift', 'Fig11_butler.png')
    figure(data, 'rel', 'Fig12_butler_equilibrium_start.png')


if __name__ == '__main__':
    main()
