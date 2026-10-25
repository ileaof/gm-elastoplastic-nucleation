"""
Surface energy of the nucleus up to and beyond 10^5 K/m, elastic and elasto-plastic Gurtin-Murdoch surface, for the
three alloys (undercooling limit of the sweep raised to GM_DT_MAX = 30 K). Marks the yield point (end of the elastic
range) and the values at 10^5 K/m.
Figure: nucleacao/figuras/energia_superficie_1e5.png.  Usage:  python energia_superficie_alto_gradiente.py
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
    os.environ.update(GM_SURFACE_MODEL=model, GM_DT_MAX='30')
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
    out = {k: np.asarray(ns[k][1:n], float) for k in ('dTdr_v', 'gam_hom_v', 'DT_v', 'rc_v')}
    s = ns['gm_surf']
    out['G_yield'] = float(out['dTdr_v'][int(np.argmin(np.abs(out['rc_v'] - s.yielded_at)))]) if s.yielded_at else None
    out['gamma_0'] = ns['gamma_0']
    return out


fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
for ax, (alloy, script, G_mean) in zip(axs, MODELS):
    e, p = run(script, 'elastic'), run(script, 'elastoplastic')
    ax.axvspan(min(e['dTdr_v'][0], p['dTdr_v'][0]), p['G_yield'], color='#2a78d6', alpha=0.08, lw=0)
    ax.plot(e['dTdr_v'], e['gam_hom_v'], color='#2a78d6', lw=2, label='Elastic GM surface')
    ax.plot(p['dTdr_v'], p['gam_hom_v'], color='#eb6834', lw=2, ls='--', label='Elasto-plastic GM surface')
    ax.axhline(e['gamma_0'], color='0.5', lw=1, ls=':', label=r'Equilibrium $\gamma_0$')
    ax.axvline(G_mean, color='0.4', lw=1, ls=':', label='Mean experimental gradient')
    ax.axvline(p['G_yield'], color='k', lw=1, ls='-.', label=f'Yield: elastic range (shaded) ends')
    ax.axvline(G_MARK, color='#b02a2a', lw=1, ls='--')
    ge, gp = np.interp(G_MARK, e['dTdr_v'], e['gam_hom_v']), np.interp(G_MARK, p['dTdr_v'], p['gam_hom_v'])
    te, tp = np.interp(G_MARK, e['dTdr_v'], e['DT_v']), np.interp(G_MARK, p['dTdr_v'], p['DT_v'])
    re_, rp = np.interp(G_MARK, e['dTdr_v'], e['rc_v']), np.interp(G_MARK, p['dTdr_v'], p['rc_v'])
    ax.plot([G_MARK, G_MARK], [ge, gp], 'o', color='#b02a2a', ms=6, zorder=5)
    ax.annotate(f'$10^5$ K/m:\nelastic {ge:.2f} J/m$^2$\nelasto-plastic {gp:.2f} J/m$^2$\n$\\Delta T$ = {tp:.2f} K, $r_C$ = {rp*1e6:.2f} $\\mu$m',
                xy=(G_MARK, gp), xytext=(0.03, 0.62), textcoords='axes fraction', fontsize=9,
                bbox=dict(boxstyle='round', fc='white', ec='0.6'), arrowprops=dict(arrowstyle='->', color='#b02a2a'))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'Thermal gradient $\nabla T$ [K m$^{-1}$]')
    ax.set_title(alloy)
    ax.grid(True, which='both', ls='--', color='0.85', lw=0.5)
    print(f'{alloy}: at 1e5 K/m gamma elastic {ge:.4f} | elasto-plastic {gp:.4f} J/m2 ({100*(gp/ge-1):+.2f} %), '
          f'gamma/gamma0 = {gp/e["gamma_0"]:.1f}, DT {te:.3f} | {tp:.3f} K, r_C {re_*1e6:.3f} | {rp*1e6:.3f} um; '
          f'sweep up to {e["dTdr_v"][-1]:.3g} | {p["dTdr_v"][-1]:.3g} K/m, yield at {p["G_yield"]:.0f} K/m')
axs[0].set_ylabel(r'Surface energy $\gamma_{SL}$ [J m$^{-2}$]')
axs[2].legend(fontsize=8.5, loc='lower right', frameon=True, edgecolor='k')
fig.tight_layout()
fig.savefig(os.path.join(FIG, 'energia_superficie_1e5.png'), dpi=200)
