"""
Curve-by-curve comparison of the elastic and of the elasto-plastic surface against the thermal gradient
dTdr_v = Delta T / (8 pi r) of the sweep: each copy runs twice (GM_SURFACE_MODEL=elastic and the default
elasto-plastic mode) and the homogeneous quantities of the script are overlaid (Gibbs-Thomson: 2nd order).
Figures: nucleacao/figuras/curvas_EP_<model>.png.  Usage:  python comparar_curvas_EP.py
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
    # alloy: (script, mean experimental gradient [K/m])
    'Al-3Cu-5Nb-0.1Fe': ('Nucleation_model_FCC_A1_new_paper_B_MRS_Meeting_2026_dev_full_GLH_AlCuNbFe_signal_beta_2nd_Sim_2026_paper_EP.py', 3294.0),
    'Al-0.8Si-0.6Mg-0.2Fe': ('Nucleation_model_FCC_A1_new_paper_2024_dev_full_GLH_AlSiMgFe_signal_beta_2nd_Sim_Luane_EP.py', 7902.38),
    'Al-3Cu-0.5Mg-0.2Fe': ('Nucleation_model_FCC_A1_new_paper_MR_2025_dev_full_GLH_AlCuMgFe_signal_beta_2nd_Sim_2026_paper_EP.py', 3708.2),
}

# (array of the script, y label)
PANELS = [
    ('tau_total', r'Surface tension $\tau_{SL}$, elastic + plastic [N/m]'),
    ('gam_hom_v', r'Surface energy $\gamma_{SL}$ [J/m$^2$]'),
    ('dfgamdr_ana_v', r'$\partial\gamma_{SL}/\partial r$ [J/m$^3$]'),
    ('GT_hom_2nd_v', r'Gibbs-Thomson 2nd order $\Gamma$ [m K]'),
    ('DT_v', r'Undercooling $\Delta T$ [K]'),
    ('rc_v', r'Critical radius $r_C$ [m]'),
    ('DGc_hom_2nd_expr_v', r'Critical free energy $\Delta G_C$ [J]'),
    ('stat_hom_2nd_v', r'Homogeneous rate $I$ [m$^{-3}$s$^{-1}$]'),
    ('stat_het_2nd_v', r'Heterogeneous rate $I$ [m$^{-3}$s$^{-1}$]'),
]


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
    out = {k: np.asarray(ns[k][1:n], dtype=float) for k in ['dTdr_v'] + [p for p, _ in PANELS] if k in ns}
    out['tau_total'] = np.asarray(ns['gm_surf'].hist[1:n], dtype=float)[:, 3]   # elastic + plastic components
    out['G_yield'] = yield_gradient(ns)
    return out, ns['gm_surf'].describe()


def yield_gradient(ns):
    """Thermal gradient at which the interface yields (end of the elastic range), None for the elastic surface."""
    ry = ns['gm_surf'].yielded_at
    if ry is None:
        return None
    n = ns['npoints']
    k = int(np.argmin(np.abs(np.asarray(ns['rc_v'][1:n]) - ry)))
    return float(ns['dTdr_v'][1:n][k])


for alloy, (script, G_mean) in MODELS.items():
    E, _ = run(script, 'elastic')
    P, desc = run(script, 'elastoplastic')
    fig, axs = plt.subplots(3, 3, figsize=(15, 12))
    for ax, (key, ylabel) in zip(axs.flat, PANELS):
        ax.plot(E['dTdr_v'], E[key], color='#2a78d6', lw=2.0)
        ax.plot(P['dTdr_v'], P[key], color='#eb6834', lw=2.0, ls='--')
        ax.axvline(G_mean, color='0.4', lw=1.0, ls=':')
        ax.axvspan(min(E['dTdr_v'][0], P['dTdr_v'][0]), P['G_yield'], color='#2a78d6', alpha=0.08, lw=0)   # elastic range
        ax.axvline(P['G_yield'], color='k', lw=1.0, ls='-.')
        ax.set_xscale('log')
        ax.set_xlabel(r'Thermal gradient $\nabla T$ [K/m]', fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.grid(True, linestyle='--', color='darkgray', linewidth=0.5)
    fig.legend(axs.flat[0].get_lines(), ['Elastic Gurtin-Murdoch surface',
                'Elasto-plastic surface, permanent deformation stored',
                f'Mean experimental gradient ({G_mean:g} K/m)',
                f'Elastic range (shaded) ends, plastic begins: yield at {P["G_yield"]:.0f} K/m'],
               loc='lower center', ncol=2, fontsize=11, frameon=True, edgecolor='k')
    fig.suptitle(f'{alloy}: elastic vs elasto-plastic Gurtin-Murdoch surface (homogeneous, 2nd order)', fontsize=14)
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    tag = alloy.replace('.', '').replace('-', '')
    fig.savefig(os.path.join(FIG, f'curvas_EP_{tag}.png'), dpi=150)
    plt.close(fig)
    k = lambda D: int(np.argmin(np.abs(D['dTdr_v'] - G_mean)))
    ke, kp = k(E), k(P)
    print(f'{alloy} at the mean gradient {G_mean:g} K/m (elastic | elasto-plastic): ' +
          ', '.join(f'{key} {E[key][ke]:.4g} | {P[key][kp]:.4g}' for key, _ in PANELS))
    print(f'   {desc}')
