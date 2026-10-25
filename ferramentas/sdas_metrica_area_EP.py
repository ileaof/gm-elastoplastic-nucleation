"""
MRB SDAS against the local solidification time, three alloys: RB and 1st order of the SDAS script (dotted line),
2nd order with the script form of the surface energy, and 2nd order with the metric of the sphere and the stress
integral over the area (initial area A0 to final area A1, excess stress s, stored energy), elastic and elasto-plastic
(stored_total, H = 0: minus the dissipated heat sigma_y alpha, alpha = |E_A| - sigma_y/|sigma_0| after yield).
Figure: sdas/figuras/sdas_metrica_area_EP.png.  Usage:  python sdas_metrica_area_EP.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, 'formas_energia_gt_sdas.py'), encoding='utf-8').read()
     .split('\nfor alloy, (script, g0, s0, dSv, G_mean) in ALLOYS.items():')[0])
SIGMA_Y = {'Al-3Cu-5Nb-0.1Fe': 0.038572, 'Al-0.8Si-0.6Mg-0.2Fe': 0.157175, 'Al-3Cu-0.5Mg-0.2Fe': 0.0422712}

fig, axs = plt.subplots(1, 3, figsize=(16, 5.4))
for ax, (alloy, (script, g0, s0, dSv, G_mean)) in zip(axs, ALLOYS.items()):
    s = 2 * s0 * np.log(y / y0)
    sy = SIGMA_Y[alloy]
    alpha = np.maximum(np.abs(s) / s0 - sy / s0, 0.0)
    area = -cum(np.abs(s) * 8 * np.pi * r) / (4 * np.pi * R0 ** 2)      # (1/A0) int_{A1}^{A0} |s| dA
    forms = [(r'2nd order, script form: $\gamma_0(r_{eq}/r)^2-s$', g0 / x - s, '#2a78d6', '-'),
             (r'2nd order, metric + area, elastic', 3 * g0 / (1 + 2 * x) + area, '#1baf7a', '-.'),
             (r'2nd order, metric + area, elasto-plastic', 3 * g0 / (1 + 2 * x) + area - sy * alpha, '#eb6834', '--')]
    ns = sdas_namespace(script)
    P, VL, tSL, xsc, ysc = ns['P'], ns['VL'], ns['tSL'], ns['xsc'], ns['ysc']
    k = (xsc > 0) & (ysc > 0)
    d = ns['diff_length_scale'](P, VL, tSL)

    def sdas(GT):
        S, *_ = ns['calc'](GT * np.ones(len(P)), *d)
        S = S * 1e6
        z = np.polyfit(np.log(tSL), np.log(S), 1)
        return S, np.exp(z[1]), z[0], 100 * np.mean(np.abs(np.interp(xsc[k], tSL, S) - ysc[k]) / ysc[k])

    src = open(os.path.join(SRC, script), encoding='utf-8').read()
    SRB = np.asarray(eval(src.split('\nSDAS_RB = ', 1)[1].split('\n', 1)[0], ns), dtype=float) * 1e6
    zr = np.polyfit(np.log(tSL), np.log(SRB), 1)
    ax.plot(tSL, SRB, color='k', lw=1.2, label=r'RB (script), SDAS = %.2f $t_{SL}^{%.3f}$' % (np.exp(zr[1]), zr[0]))
    S1, c1, e1, m1 = sdas(ns['Gibbs_Thomson_het_1st'])
    ax.plot(tSL, S1, color='k', ls=':', lw=2.2,
            label=r'1st order (script), $\Gamma$ = %.3e: %.2f $t_{SL}^{%.3f}$, MAPE %.1f %%' % (ns['Gibbs_Thomson_het_1st'], c1, e1, m1))
    for lab, g, col, ls in forms:
        G = g / (4 * np.pi * r * r * dSv)
        m = np.r_[True, np.diff(G) > 0].cumprod().astype(bool)
        G2 = np.interp(G_mean, G[m], g[m]) / dSv
        S, c, e, mp = sdas(G2)
        ax.plot(tSL, S, color=col, ls=ls, lw=2,
                label=lab + '\n' + r'$\Gamma^{2nd}$ = %.3e: %.2f $t_{SL}^{%.3f}$, MAPE %.1f %%' % (G2, c, e, mp))
    yerr = [np.abs(ns['yscm'][k]), np.abs(ns['yscp'][k])] if ('yscm' in ns and 'yscp' in ns) else None
    ax.errorbar(xsc[k], ysc[k], yerr=yerr, fmt='sk', ms=5, capsize=2, lw=1, label='Experimental')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'Local solidification time, $t_{SL}$ [s]')
    ax.set_title(alloy + r' ($\sigma_y^s$ = %g N/m)' % sy, fontsize=11)
    ax.grid(True, which='both', ls='--', color='0.85', lw=0.5)
    ax.legend(fontsize=6.8, loc='upper left', frameon=True, edgecolor='k')
axs[0].set_ylabel(r'Secondary dendrite arm spacing, SDAS [$\mu$m]')
fig.tight_layout()
out = os.path.join(HERE, '..', 'sdas', 'figuras', 'sdas_metrica_area_EP.png')
fig.savefig(out, dpi=200)
print('figure', os.path.abspath(out))
