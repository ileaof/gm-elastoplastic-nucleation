"""
Surface tension composed of the elastic and plastic components, tau = (sigma_0 + s^e) + |sigma_0| E_A^p + sigma_y
(H = 0), up to and beyond 10^5 K/m, against the elastic Gurtin-Murdoch surface tension, for the three alloys.
Marks the yield point (end of the elastic range) and the values at 10^5 K/m.
Figure: nucleacao/figuras/tensao_superficial_1e5.png.  Usage:  python tensao_superficial_alto_gradiente.py
"""
import os
import sys
import runpy
import tempfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
NUC = os.path.join(HERE, '..', 'nucleacao')
FIG = os.path.join(NUC, 'figuras')
sys.path.insert(0, NUC)
G_MARK = 1e5

MODELS = [
    ('Al-3Cu-5Nb-0.1Fe', 'Nucleation_model_FCC_A1_new_paper_B_MRS_Meeting_2026_dev_full_GLH_AlCuNbFe_signal_beta_2nd_Sim_2026_paper_EP.py', 3294.0),
    ('Al-0.8Si-0.6Mg-0.2Fe', 'Nucleation_model_FCC_A1_new_paper_2024_dev_full_GLH_AlSiMgFe_signal_beta_2nd_Sim_Luane_EP.py', 7902.38),
    ('Al-3Cu-0.5Mg-0.2Fe', 'Nucleation_model_FCC_A1_new_paper_MR_2025_dev_full_GLH_AlCuMgFe_signal_beta_2nd_Sim_2026_paper_EP.py', 3708.2),
]


def run(script, model):
    os.environ.update(GM_SURFACE_MODEL=model, GM_DT_MAX='30', GM_DUMP_RATE='')
    plt.show = lambda *a, **k: plt.close('all')
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as wd:
        os.chdir(wd)
        try:
            ns = runpy.run_path(os.path.join(NUC, script), run_name='__main__')
        finally:
            os.chdir(cwd)
            plt.close('all')
    n = ns['npoints']
    h = np.asarray(ns['gm_surf'].hist[1:n], float)
    out = {'G': np.asarray(ns['dTdr_v'][1:n], float), 'rc': np.asarray(ns['rc_v'][1:n], float),
           'tau_e': h[:, 0], 'tau_p': h[:, 1], 'tau': h[:, 3]}
    s = ns['gm_surf']
    out['G_yield'] = float(out['G'][int(np.argmin(np.abs(out['rc'] - s.yielded_at)))]) if s.yielded_at else None
    return out


fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
for ax, (alloy, script, G_mean) in zip(axs, MODELS):
    e, p = run(script, 'elastic'), run(script, 'elastoplastic')
    ax.axvspan(min(e['G'][0], p['G'][0]), p['G_yield'], color='#2a78d6', alpha=0.08, lw=0)
    ax.plot(e['G'], e['tau'], color='#2a78d6', lw=2, label='Elastic GM surface')
    ax.plot(p['G'], p['tau'], color='#eb6834', lw=2, ls='--', zorder=3,
            label=r'Elasto-plastic, total: $\sigma_0+s^e+|\sigma_0|E_A^p+\sigma_y^s$')
    ax.plot(p['G'], p['tau_e'], color='k', lw=1.2, ls=':', label=r'Elastic component, $\sigma_0+s^e$')
    ax.plot(p['G'], p['tau_p'], color='#8a3ffc', lw=1.2, ls='-.', label=r'Plastic component, $|\sigma_0|E_A^p$')
    ax.axvline(G_mean, color='0.4', lw=1, ls=':')
    ax.axvline(p['G_yield'], color='k', lw=1, ls='-.', label='Yield: elastic range (shaded) ends')
    ax.axvline(G_MARK, color='#b02a2a', lw=1, ls='--')
    vals = {}
    for tag, G in (('mean', G_mean), ('1e5', G_MARK)):
        vals[tag] = [np.interp(G, e['G'], e['tau']), np.interp(G, p['G'], p['tau']),
                     np.interp(G, p['G'], p['tau_e']), np.interp(G, p['G'], p['tau_p'])]
    te, tp = vals['1e5'][0], vals['1e5'][1]
    ax.plot([G_MARK, G_MARK], [te, tp], 'o', color='#b02a2a', ms=6, zorder=5)
    ax.text(0.40, 0.97, f'$10^5$ K/m:\nelastic {te:.3f} N/m\nelasto-plastic {tp:.3f} N/m', transform=ax.transAxes,
            va='top', fontsize=9, bbox=dict(boxstyle='round', fc='white', ec='0.6'))
    ax.set_xscale('log')
    ax.set_xlabel(r'Thermal gradient $\nabla T$ [K m$^{-1}$]')
    ax.set_title(alloy)
    ax.grid(True, which='both', ls='--', color='0.85', lw=0.5)
    for tag in ('mean', '1e5'):
        v = vals[tag]
        print(f'{alloy} at {tag}: tau elastic {v[0]:.4f} | EP total {v[1]:.4f} (elastic comp {v[2]:.4f}, '
              f'plastic comp {v[3]:.4f}) N/m, {100 * (v[1] / v[0] - 1):+.2f} %')
axs[0].set_ylabel(r'Surface tension $\tau_{SL}$ [N m$^{-1}$]')
axs[0].legend(fontsize=8, loc='lower left', frameon=True, edgecolor='k')
fig.tight_layout()
fig.savefig(os.path.join(FIG, 'tensao_superficial_1e5.png'), dpi=200)
