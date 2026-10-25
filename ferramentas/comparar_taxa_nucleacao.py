"""
Nucleation rate of the elastic and of the elasto-plastic surface, same script and same Delta T range:
each copy runs twice (GM_SURFACE_MODEL=elastic and the default elasto-plastic mode) and the
homogeneous 2nd-order statistics stat_hom_2nd_v of the script are overlaid against the thermal
gradient dTdr_v = Delta T / (8 pi r) of the sweep.
Figures: nucleacao/figuras/taxa_nucleacao_<model>.png.  Usage:  python comparar_taxa_nucleacao.py
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
os.makedirs(FIG, exist_ok=True)
sys.path.insert(0, NUC)

MODELS = {
    'Al-3Cu-5Nb-0.1Fe': 'Nucleation_model_FCC_A1_new_paper_B_MRS_Meeting_2026_dev_full_GLH_AlCuNbFe_signal_beta_2nd_Sim_2026_paper_EP.py',
    'Al-0.8Si-0.6Mg-0.2Fe': 'Nucleation_model_FCC_A1_new_paper_2024_dev_full_GLH_AlSiMgFe_signal_beta_2nd_Sim_Luane_EP.py',
    'Al-3Cu-0.5Mg-0.2Fe': 'Nucleation_model_FCC_A1_new_paper_MR_2025_dev_full_GLH_AlCuMgFe_signal_beta_2nd_Sim_2026_paper_EP.py',
}


def run(script, surface_model):
    env_old = os.environ.get('GM_SURFACE_MODEL')
    os.environ['GM_SURFACE_MODEL'] = surface_model
    plt.show = lambda *a, **k: plt.close('all')
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as wd:
        os.chdir(wd)
        try:
            ns = runpy.run_path(os.path.join(NUC, script), run_name='__main__')
        finally:
            os.chdir(cwd)
            plt.close('all')
            if env_old is None:
                os.environ.pop('GM_SURFACE_MODEL', None)
            else:
                os.environ['GM_SURFACE_MODEL'] = env_old
    n = ns['npoints']
    G_y = None
    if ns['gm_surf'].yielded_at is not None:   # thermal gradient at the yield point (end of the elastic range)
        G_y = float(ns['dTdr_v'][1:n][int(np.argmin(np.abs(np.asarray(ns['rc_v'][1:n]) - ns['gm_surf'].yielded_at)))])
    return ns['dTdr_v'][1:n], ns['DT_v'][1:n], ns['stat_hom_2nd_v'][1:n], (ns['gm_surf'].describe(), G_y)


for alloy, script in MODELS.items():
    G_e, DT_e, I_e, _ = run(script, 'elastic')
    G_p, DT_p, I_p, (desc, G_y) = run(script, 'elastoplastic')
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.plot(G_e, I_e, color='#2a78d6', lw=2.0, label='Elastic Gurtin-Murdoch surface')
    ax.plot(G_p, I_p, color='#eb6834', lw=2.0, ls='--',
            label='Elasto-plastic surface, permanent deformation stored')
    ax.axvspan(min(G_e[0], G_p[0]), G_y, color='#2a78d6', alpha=0.08, lw=0)
    ax.axvline(G_y, color='k', lw=1.0, ls='-.', label=f'Elastic range (shaded) ends, plastic begins: yield at {G_y:.0f} K/m')
    ax.set_xscale('log')
    ax.set_xlabel(r'Thermal gradient $\nabla T$ [K/m]', fontsize=13)
    ax.set_ylabel(r'Nucleation Sites, $I$ [m$^{-3}$s$^{-1}$]', fontsize=13)
    ax.set_title(f'{alloy}: homogeneous 2nd-order statistics', fontsize=12)
    ax.grid(True, linestyle='--', color='darkgray', linewidth=0.5)
    ax.legend(fontsize=10, frameon=True, edgecolor='k', loc='upper center', bbox_to_anchor=(0.5, -0.15))
    tag = alloy.replace('.', '').replace('-', '')
    fig.savefig(os.path.join(FIG, f'taxa_nucleacao_{tag}.png'), dpi=200, bbox_inches='tight')
    ke, kp = int(np.argmax(I_e)), int(np.argmax(I_p))
    print(f'{alloy}: elastic peak I = {I_e[ke]:.3e} at grad T = {G_e[ke]:.4g} K/m (DT = {DT_e[ke]:.3f} K) | '
          f'elasto-plastic peak I = {I_p[kp]:.3e} at grad T = {G_p[kp]:.4g} K/m (DT = {DT_p[kp]:.3f} K) | {desc}')
