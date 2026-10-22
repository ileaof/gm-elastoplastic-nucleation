"""
Curve-by-curve comparison of the elastic and elasto-plastic Gurtin-Murdoch (GM) surfaces against the
thermal gradient, Al-3Cu-5Nb-0.1Fe (FCC_A1), from the tablefull_* files written by
Nucleation_model_..._ElastoPlastic_GM.py:

    homogeneous surface tension   tau = -sigma_0 + s      (column 'sigma')
    homogeneous surface energy    gamma_SL                (column 'gam_hom')
    Gibbs-Thomson 1st order       Gamma^1st               (column 'GT_het')
    Gibbs-Thomson 2nd order       Gamma^2nd               (column 'GT_het_2nd')

Run the nucleation model first for each case:
    GM_SURFACE_MODEL=elastic python Nucleation_model_..._ElastoPlastic_GM.py
    python Nucleation_model_..._ElastoPlastic_GM.py
    GM_PLASTIC_REFERENCE=1 GM_SIGMA_YA0=0.0695 GM_H_A=0 python Nucleation_model_..._ElastoPlastic_GM.py
"""
import csv
import numpy as np
import matplotlib.pyplot as plt

cases = [
    ('Elastic GM', 'tablefull_paper_B_MRS_2026_GM_elastic.csv', '#2a78d6', '-'),
    (r'Elasto-plastic GM, fixed reference ($\sigma_y^s$ = 0.5 N/m, $H$ = 0.1$\sigma_0$)',
     'tablefull_paper_B_MRS_2026_GM_elastoplastic_stored.csv', '#eb6834', '--'),
    (r'Elasto-plastic GM, updated reference ($\sigma_y^s$ = 0.0695 N/m, $H$ = 0)',
     'tablefull_paper_B_MRS_2026_GM_elastoplastic_stored_refupd_sy0.0695_H0.csv', '#1baf7a', '-.'),
]
quantities = [
    ('sigma',      r'Surface tension, $\tau = -\sigma_0 + s$ [N m$^{-1}$]', 'linear'),
    ('gam_hom',    r'Surface energy, $\gamma_{SL}$ [J m$^{-2}$]',            'log'),
    ('GT_het',     r'Gibbs-Thomson 1st order, $\Gamma^{1st}$ [m K]',         'log'),
    ('GT_het_2nd', r'Gibbs-Thomson 2nd order, $\Gamma^{2nd}$ [m K]',         'log'),
]

# experimental thermal gradients G_L = T_R / V_L (Thermal_data_Ivaldo_B_MRS_Meeting_2026_Al3Cu5Nb_GL_data_indev.py)
POS = np.array([5.0, 10.0, 14.0, 19.0, 29.0, 50.0, 69.0, 90.0])
TR = np.array([11.81538, 7.58399, 7.67473, 4.65463, 2.53968, 1.04679, 0.33553, 0.29251])
VL = np.array([1.31773, 0.89396, 0.71244, 0.60647, 0.48333, 0.36312, 0.30079, 0.26131])
G_EXP = TR / VL * 1000.0
G_MEAN = 3294.0
G_MAX_PLOT = 1.1e5


def load(fname):
    rows = list(csv.reader(open(fname), delimiter=';'))
    head = [h.strip() for h in rows[0]]
    data = [r for r in rows[1:] if not r[0].strip().endswith('.5')]      # drop the interpolated row
    col = lambda name: np.array([float(r[head.index(name)]) for r in data])
    out = {q: col(q) for q, _, _ in quantities}
    out['G'] = col('dTdr')
    out['DT'] = col('DT')
    out['r'] = col('r_hom_2nd')
    return out


def at(d, q, G):
    ok = np.isfinite(d[q])
    if not ok.any() or G < d['G'][ok][0] or G > d['G'][ok][-1]:
        return float('nan')
    return float(np.interp(G, d['G'][ok], d[q][ok]))


data = [load(f) for _, f, _, _ in cases]

# table at the experimental gradients and at the mean gradient
tags = ['elastic', 'EP_fixed_ref', 'EP_updated_ref']
with open('GM_elastic_vs_elastoplastic_curves_B_MRS_2026.csv', 'w') as fh:
    fh.write('P_mm;GL;' + ';'.join(f'{q}_{t}' for q, _, _ in quantities for t in tags) + '\n')
    for P_, G_ in list(zip(POS, G_EXP)) + [(float('nan'), G_MEAN)]:
        fh.write(f'{P_:g};{G_:.1f};' + ';'.join(f'{at(d, q, G_):.6g}' for q, _, _ in quantities for d in data) + '\n')

print(f'{"":>10}' + ''.join(f'{t:>16}' for t in tags) + f'{"EP_fix/el":>11}{"EP_upd/el":>11}')
for q, lab, _ in quantities:
    print(f'--- {q}')
    for P_, G_ in list(zip(POS, G_EXP)) + [(None, G_MEAN)]:
        v = [at(d, q, G_) for d in data]
        name = f'G={G_:.0f}' if P_ is not None else f'mean {G_:.0f}'
        print(f'{name:>10}' + ''.join(f'{x:>16.5g}' for x in v) + f'{v[1]/v[0]:>11.3f}{v[2]/v[0]:>11.3f}')

fig, axs = plt.subplots(2, 2, figsize=(12.5, 9.0), sharex=True)
for ax, (q, lab, yscale) in zip(axs.flat, quantities):
    for (clab, _, col, ls), d in zip(cases, data):
        m = np.isfinite(d[q]) & (d['G'] <= G_MAX_PLOT)
        ax.plot(d['G'][m], d[q][m], color=col, linestyle=ls, linewidth=2.0, label=clab)
    ax.axvspan(G_EXP.min(), G_EXP.max(), color='0.85', alpha=0.5, lw=0,
               label=r'Experimental $\nabla T_L$ range (%.0f-%.0f K m$^{-1}$)' % (G_EXP.min(), G_EXP.max()))
    ax.axvline(G_MEAN, color='k', linestyle=':', lw=1.2, label=r'Mean $\overline{\nabla T}$ = %g K m$^{-1}$' % G_MEAN)
    ax.set_xscale('log')
    ax.set_yscale(yscale)
    ax.set_ylabel(lab, fontsize=12)
    ax.grid(True, which='both', linestyle='--', color='darkgray', linewidth=0.5)
for ax in axs[1]:
    ax.set_xlabel(r'Thermal gradient, $\nabla T$ [K m$^{-1}$]', fontsize=12)
axs[0, 0].set_title('(a) homogeneous surface tension', fontsize=12, loc='left')
axs[0, 1].set_title('(b) homogeneous surface energy', fontsize=12, loc='left')
axs[1, 0].set_title('(c) Gibbs-Thomson, 1st order', fontsize=12, loc='left')
axs[1, 1].set_title('(d) Gibbs-Thomson, 2nd order', fontsize=12, loc='left')
handles, labels = axs[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=2, frameon=True, edgecolor='k', fontsize=10, bbox_to_anchor=(0.5, -0.12))
fig.suptitle('Elastic x elasto-plastic Gurtin-Murdoch surface, Al-3Cu-5Nb-0.1Fe (FCC_A1)', fontsize=14)
fig.tight_layout()
fig.savefig('Figure_GM_elastic_vs_elastoplastic_curves_B_MRS_2026.png', dpi=200, bbox_inches='tight')
