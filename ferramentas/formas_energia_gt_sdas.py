"""
Effect of the form of the surface free energy (metric of the sphere, Materials paper Eqs. 89, 101, 102) on the
1st- and 2nd-order Gibbs-Thomson tensors at the mean experimental gradient and on the MRB SDAS of the three alloys.
Elastic Gurtin-Murdoch surface (spherical solution of the scripts for the excess stress s). Only reads the originals.

forms of gamma(r), x = r^2/r0^2:
  scripts   : gamma_0/x - s                                   (before the generalization)
  eq102     : gamma_0 * 3/(1 + 2x) - s                        (Eq. 89 with the potential of Eq. 102, theta = theta_0)
  eq101     : gamma_0 * 3/(1 + 2x) - int_{r0}^{r} x ds        (free energy integrated from the potential of Eq. 101)
  dA        : gamma_0 * 3/(1 + 2x) - (1/4 pi r0^2) int s dA   (slide 13 read with the area as variable, dA = 8 pi r dr)
Gamma_2nd = gamma/|DS_V| at the gradient grad T = gamma/(4 pi r^2 |DS_V|);
Gamma_1st = Gamma_2nd / (1 + r |d gamma/dr| / (2 gamma)).
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SRC = os.path.join(ROOT, 'Sim Luane Final')
sys.path.insert(0, SRC)

LAM, MU, L0, M0, R0 = 50.34e9, 25.94e9, -5.0, 1.6, 6.3e-6
ALLOYS = {
    # alloy: (SDAS script, gamma_0, |sigma_0|, |DS_V|, mean gradient)
    'Al-3Cu-5Nb-0.1Fe': ('scriptsdasmodelAl3Cu5Nb01Fe_zoom_horizontal_B_MRS_Meeting_2026.py', 0.296, 1.3761, 979873.06, 3294.0),
    'Al-0.8Si-0.6Mg-0.2Fe': ('scriptsdasmodelAl08Si06Mg02Fe_zoom_horizontal_new_script_Ferreira_2024.py', 0.183, 1.09, 926093.58, 7902.38),
    'Al-3Cu-0.5Mg-0.2Fe': ('scriptsdasmodelAl3Cu05Mg02Fe_zoom_horizontal_new_script_Ferreira_2025.py', 0.169, 1.209, 834686.60, 3708.2),
}


def sdas_namespace(script):
    src = open(os.path.join(SRC, script), encoding='utf-8').read()
    ns = {}
    cwd = os.getcwd()
    os.chdir(SRC)
    try:
        exec(compile(src.split('\nSDAS_RB = calc_RB')[0], script, 'exec'), ns)
    finally:
        os.chdir(cwd)
    return ns


r = np.linspace(R0, 0.15 * R0, 40000)
x = (r / R0) ** 2
y = (3 * LAM + 2 * MU) * r + 2 * (L0 + M0)
y0 = (3 * LAM + 2 * MU) * R0 + 2 * (L0 + M0)
cum = lambda f: np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(r))])

for alloy, (script, g0, s0, dSv, G_mean) in ALLOYS.items():
    s = 2 * s0 * np.log(y / y0)
    ds = np.gradient(s, r)
    forms = {'scripts': g0 / x - s,
             'eq102': 3 * g0 / (1 + 2 * x) - s,
             'eq101': 3 * g0 / (1 + 2 * x) - cum(x * ds),
             'dA': 3 * g0 / (1 + 2 * x) - cum(s * 2 * r / R0 ** 2)}
    ns = sdas_namespace(script)
    P, VL, tSL, xsc, ysc = ns['P'], ns['VL'], ns['tSL'], ns['xsc'], ns['ysc']
    k = (xsc > 0) & (ysc > 0)
    d = ns['diff_length_scale'](P, VL, tSL)
    print(f'== {alloy} (mean gradient {G_mean:g} K/m)')
    for name, g in forms.items():
        G = g / (4 * np.pi * r * r * dSv)
        m = np.r_[True, np.diff(G) > 0].cumprod().astype(bool)
        if G[m].max() < G_mean:
            print(f'   {name:8s} does not reach the mean gradient (max {G[m].max():.0f} K/m)')
            continue
        dg = np.gradient(g, r)
        gm = np.interp(G_mean, G[m], g[m])
        rm = np.interp(G_mean, G[m], r[m])
        dgm = np.interp(G_mean, G[m], dg[m])
        G2 = gm / dSv
        G1 = G2 / (1 + rm * abs(dgm) / (2 * gm))
        out = []
        for GT in (G2, G1):
            S, *_ = ns['calc'](GT * np.ones(len(P)), *d)
            S = S * 1e6
            z = np.polyfit(np.log(tSL), np.log(S), 1)
            Se = np.interp(xsc[k], tSL, S)
            out.append((np.exp(z[1]), z[0], 100 * np.mean(np.abs(Se - ysc[k]) / ysc[k]),
                        np.mean(np.log(Se / ysc[k]))))
        print(f'   {name:8s} gamma {gm:.4f} J/m2 | Gamma_2nd {G2:.4e} m K: SDAS {out[0][0]:.2f} t^{out[0][1]:.3f}, '
              f'MAPE {out[0][2]:.1f} %, bias {out[0][3]:+.3f} | Gamma_1st {G1:.4e} m K: SDAS {out[1][0]:.2f} '
              f't^{out[1][1]:.3f}, MAPE {out[1][2]:.1f} %, bias {out[1][3]:+.3f}')


# ---------------------------------------------------------------- SDAS against t_SL, 2nd order, each form
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

STYLE = {'scripts': (r'Script form (traditional): $\gamma_0(r_{eq}/r)^2-s$', '#2a78d6', '-'),
         'eq102': (r'New paper, Eqs. (89)+(102): $3\gamma_0/(1+2r^2/r_0^2)-s$', '#eb6834', '--'),
         'eq101': (r'New paper, Eqs. (89)+(101): $3\gamma_0/(1+2r^2/r_0^2)-\int (r^2/r_0^2)\,ds$', '#1baf7a', '-.')}
fig, axs = plt.subplots(1, 3, figsize=(16, 5.2))
for ax, (alloy, (script, g0, s0, dSv, G_mean)) in zip(axs, ALLOYS.items()):
    s = 2 * s0 * np.log(y / y0)
    ds = np.gradient(s, r)
    forms = {'scripts': g0 / x - s, 'eq102': 3 * g0 / (1 + 2 * x) - s, 'eq101': 3 * g0 / (1 + 2 * x) - cum(x * ds)}
    ns = sdas_namespace(script)
    P, VL, tSL, xsc, ysc = ns['P'], ns['VL'], ns['tSL'], ns['xsc'], ns['ysc']
    k = (xsc > 0) & (ysc > 0)
    d = ns['diff_length_scale'](P, VL, tSL)
    for name, g in forms.items():
        G = g / (4 * np.pi * r * r * dSv)
        m = np.r_[True, np.diff(G) > 0].cumprod().astype(bool)
        G2 = np.interp(G_mean, G[m], g[m]) / dSv
        S, *_ = ns['calc'](G2 * np.ones(len(P)), *d)
        S = S * 1e6
        z = np.polyfit(np.log(tSL), np.log(S), 1)
        lab, col, ls = STYLE[name]
        ax.plot(tSL, S, color=col, ls=ls, lw=2,
                label=lab + '\n' + r'$\Gamma^{2nd}$ = %.3e m K, SDAS = %.2f $t_{SL}^{%.3f}$' % (G2, np.exp(z[1]), z[0]))
    # the normal curves of the SDAS script: MRB with its 1st-order Gamma and the RB model (script constants)
    S1, *_ = ns['calc'](ns['Gibbs_Thomson_het_1st'] * np.ones(len(P)), *d)
    S1 = S1 * 1e6
    z1 = np.polyfit(np.log(tSL), np.log(S1), 1)
    ax.plot(tSL, S1, color='k', ls=':', lw=2,
            label='Script: MRB, $\\Gamma^{1st}$ = %.3e m K, SDAS = %.2f $t_{SL}^{%.3f}$'
                  % (ns['Gibbs_Thomson_het_1st'], np.exp(z1[1]), z1[0]))
    src = open(os.path.join(SRC, script), encoding='utf-8').read()
    call = src.split('\nSDAS_RB = ', 1)[1].split('\n', 1)[0]
    SRB = np.asarray(eval(call, ns), dtype=float) * 1e6
    zr = np.polyfit(np.log(tSL), np.log(SRB), 1)
    ax.plot(tSL, SRB, color='k', ls='-', lw=1.2,
            label='Script: RB (Rappaz and Boettinger), SDAS = %.2f $t_{SL}^{%.3f}$' % (np.exp(zr[1]), zr[0]))
    yerr = None
    if 'yscm' in ns and 'yscp' in ns:
        yerr = [np.abs(ns['yscm'][k]), np.abs(ns['yscp'][k])]
    ax.errorbar(xsc[k], ysc[k], yerr=yerr, fmt='sk', ms=5, capsize=2, lw=1, label='Experimental')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'Local solidification time, $t_{SL}$ [s]')
    ax.set_title(alloy + r' (MRB, 2nd order, mean $\nabla T$ = %g K/m)' % G_mean, fontsize=10)
    ax.grid(True, which='both', ls='--', color='0.85', lw=0.5)
    ax.legend(fontsize=7, loc='upper left', frameon=True, edgecolor='k')
axs[0].set_ylabel(r'Secondary dendrite arm spacing, SDAS [$\mu$m]')
fig.tight_layout()
out = os.path.join(HERE, '..', 'sdas', 'figuras', 'sdas_formas_energia_2a_ordem.png')
fig.savefig(out, dpi=200)
print('figure', out)
