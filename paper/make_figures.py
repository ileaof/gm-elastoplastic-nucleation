"""Figures of the elasto-plastic surface paper (rate of the three alloys, restoration along the casting)."""
import os, sys, runpy, tempfile
import numpy as np
from scipy import integrate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

EP = r'C:\MacbookPro\ufpa\PPGEM\SDAS Al_Si_Mg_Fe\modelos_elasto-plastos'
NUC = os.path.join(EP, 'nucleacao')
OUT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, NUC)
plt.rcParams.update({'font.size': 11, 'font.family': 'serif', 'mathtext.fontset': 'dejavuserif'})

MODELS = [
    ('Al-3Cu-5Nb-0.1Fe', 'Nucleation_model_FCC_A1_new_paper_B_MRS_Meeting_2026_dev_full_GLH_AlCuNbFe_signal_beta_2nd_Sim_2026_paper_EP.py', 3294.0),
    ('Al-0.8Si-0.6Mg-0.2Fe', 'Nucleation_model_FCC_A1_new_paper_2024_dev_full_GLH_AlSiMgFe_signal_beta_2nd_Sim_Luane_EP.py', 7902.38),
    ('Al-3Cu-0.5Mg-0.2Fe', 'Nucleation_model_FCC_A1_new_paper_MR_2025_dev_full_GLH_AlCuMgFe_signal_beta_2nd_Sim_2026_paper_EP.py', 3708.2),
]


def run(script, model):
    os.environ['GM_SURFACE_MODEL'] = model
    show = plt.show
    plt.show = lambda *a, **k: plt.close('all')
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as wd:
        os.chdir(wd)
        try:
            ns = runpy.run_path(os.path.join(NUC, script), run_name='__main__')
        finally:
            os.chdir(cwd)
            plt.close('all')
            plt.show = show
    n = ns['npoints']
    out = {k: np.asarray(ns[k][1:n], float) for k in ('dTdr_v', 'stat_hom_2nd_v', 'GT_het_2nd_v', 'rc_v')}
    s = ns['gm_surf']
    out['G_yield'] = float(out['dTdr_v'][int(np.argmin(np.abs(out['rc_v'] - s.yielded_at)))]) if s.yielded_at else None
    out['sigma_y'], out['gamma_0'], out['DSv'] = s.sigma_y, ns['gamma_0'], abs(ns['Delta_Sv_hom'])
    return out


fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.2))
runs = {}
for ax, (alloy, script, G) in zip(axs, MODELS):
    e, p = run(script, 'elastic'), run(script, 'elastoplastic')
    runs[alloy] = (e, p)
    ax.axvspan(min(e['dTdr_v'][0], p['dTdr_v'][0]), p['G_yield'], color='#2a78d6', alpha=0.08, lw=0)
    ax.plot(e['dTdr_v'], e['stat_hom_2nd_v'], color='#2a78d6', lw=2, label='Elastic GM surface')
    ax.plot(p['dTdr_v'], p['stat_hom_2nd_v'], color='#eb6834', lw=2, ls='--', label='Elasto-plastic GM surface')
    ax.axvline(G, color='0.4', lw=1, ls=':', label=r'Mean experimental $\nabla T$')
    ax.axvline(p['G_yield'], color='k', lw=1, ls='-.', label='Yield: elastic range (shaded) ends')
    ax.set_xscale('log')
    ax.set_xlabel(r'Thermal gradient $\nabla T$ [K m$^{-1}$]')
    ax.set_title(alloy, fontsize=11)
    ax.grid(True, ls='--', color='0.8', lw=0.5)
axs[0].set_ylabel(r'Nucleation rate $I$ [m$^{-3}$ s$^{-1}$]')
axs[2].legend(fontsize=8.5, loc='upper right', frameon=True, edgecolor='k')
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'Fig2_nucleation_rate.png'), dpi=300)
plt.close(fig)

# elastic restoration along the casting of Al-3Cu-5Nb-0.1Fe: only the recoverable part |s| <= sigma_y is restored
_, p = runs['Al-3Cu-5Nb-0.1Fe']
P = np.linspace(5, 90, 171)
GL = 10405.54 - 200.22 * P + 1.04692 * P ** 2              # quadratic fit of the thermal data (5-90 mm)
VL = 3.4 * P ** -0.56 * 1e-3                                # [m/s]
tSL = 4.93116 * P ** 0.87606
Tdot = GL * VL * 1e3                                        # GL [K/m] * VL [mm/s] * 1e-3 -> K/s (VL above already in m/s)
Tdot = GL * 3.4 * P ** -0.56 * 1e-3
DT_SL = integrate.simpson(Tdot * tSL, x=P) / (90 - 5)
beta_s = (2737.04 - 2527.05) / 2737.04
eps = beta_s * Tdot / DT_SL
GEP = np.interp(GL, p['dTdr_v'], p['GT_het_2nd_v'])
dG_el = p['sigma_y'] / p['DSv']                             # largest restorable part of Gamma
G0 = p['gamma_0'] / p['DSv']
print(f'DT_SL = {DT_SL:.2f} K, beta_s = {beta_s:.4f}, sigma_y/|DSv| = {dG_el:.4e} m K, gamma_0/|DSv| = {G0:.4e} m K')
for Pk in (5, 20, 40, 60, 80, 90):
    k = int(np.argmin(np.abs(P - Pk)))
    print(f'P={Pk:3d} GL={GL[k]:7.0f} Tdot={Tdot[k]:6.3f} eps={eps[k]:.3e} G_EP={GEP[k]:.4e} G_R(f=1)={GEP[k]-dG_el:.4e} '
          f'ratio_eq={(GEP[k]-dG_el)/G0:.2f} red={100*dG_el/GEP[k]:.1f}%')
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(P, eps, color='#2a78d6', lw=2)
ax[0].set_yscale('log')
ax[0].set_ylabel(r'Rate of density change $\dot\varepsilon_v$ [s$^{-1}$]')
ax[0].set_xlabel('Position from the cooled surface, P [mm]')
ax[1].plot(P, GEP * 1e6, color='#eb6834', lw=2, ls='--', label=r'Elasto-plastic, loaded ($f=0$)')
ax[1].fill_between(P, (GEP - dG_el) * 1e6, GEP * 1e6, color='#eb6834', alpha=0.25, lw=0,
                   label=r'Range of the elastic restoration, $0\leq f\leq 1$')
ax[1].plot(P, (GEP - dG_el) * 1e6, color='k', lw=1.5, label=r'Elastic part fully restored ($f=1$)')
ax[1].axhline(G0 * 1e6, color='0.4', lw=1, ls=':', label=r'Equilibrium, $\gamma_0/|\Delta S_V|$')
ax[1].set_xlabel('Position from the cooled surface, P [mm]')
ax[1].set_ylabel(r'Second-order thermal-field tensor $\Gamma$ [$\mu$m K]')
ax[1].set_ylim(0, None)
ax[1].legend(fontsize=9, frameon=True, edgecolor='k')
for a in ax:
    a.grid(True, ls='--', color='0.8', lw=0.5)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'Fig3_restoration.png'), dpi=300)
print('ok')
