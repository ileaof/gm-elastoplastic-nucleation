"""Point-to-point restoration fraction needed by each SDAS point of Al-3Cu-5Nb-0.1Fe."""
import os, sys, csv, runpy, tempfile
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = r'C:\MacbookPro\ufpa\PPGEM\SDAS Al_Si_Mg_Fe'
SRC = os.path.join(ROOT, 'Sim Luane Final')
EP = os.path.join(ROOT, 'modelos_elasto-plastos')
NUC = os.path.join(EP, 'nucleacao')
EST = os.path.join(EP, 'estudo_AlCuNbFe')
sys.path.insert(0, SRC); sys.path.insert(0, NUC)

# SDAS model (original script, definitions and data only)
src = open(os.path.join(SRC, 'scriptsdasmodelAl3Cu5Nb01Fe_zoom_horizontal_B_MRS_Meeting_2026.py'), encoding='utf-8').read()
ns = {}; cwd = os.getcwd(); os.chdir(SRC)
exec(compile(src.split('\nSDAS_RB = calc_RB')[0], 'sdas', 'exec'), ns); os.chdir(cwd)
P, VL, tSL, xsc, ysc = ns['P'], ns['VL'], ns['tSL'], ns['xsc'], ns['ysc']
d = ns['diff_length_scale'](P, VL, tSL)
n_exp = 16
t_exp, s_exp = xsc[:n_exp], ysc[:n_exp]

def sdas_at(G, t):
    s, *_ = ns['calc'](G * np.ones(len(P)), *d)
    return np.interp(t, tSL, s * 1e6)

G_req = np.array([brentq(lambda G: sdas_at(G, t) - s, 1e-9, 1e-4, xtol=1e-15) for t, s in zip(t_exp, s_exp)])

# local gradient (quadratic fit 5-90 mm of the thermal data, constant outside)
Pf, GLf = [], []
for row in csv.DictReader(open(os.path.join(EST, 'GT_point_to_point_B_MRS_2026_GM_elastic.csv')), delimiter=';'):
    if row['source'] == 'fit':
        Pf.append(float(row['P_mm'])); GLf.append(float(row['GL']))
P_exp = (t_exp / 4.93116) ** (1 / 0.87606)
GL_exp = np.interp(P_exp, Pf, GLf)

# 2nd-order Gamma of the nucleation sweep, elastic and elasto-plastic (complement)
def sweep(model):
    os.environ['GM_SURFACE_MODEL'] = model
    plt.show = lambda *a, **k: plt.close('all')
    with tempfile.TemporaryDirectory() as wd:
        os.chdir(wd)
        try:
            n = runpy.run_path(os.path.join(NUC, 'Nucleation_model_FCC_A1_new_paper_B_MRS_Meeting_2026_dev_full_GLH_AlCuNbFe_signal_beta_2nd_Sim_2026_paper_EP.py'), run_name='__main__')
        finally:
            os.chdir(cwd); plt.close('all')
    k = n['npoints']
    return np.asarray(n['dTdr_v'][1:k]), np.asarray(n['GT_hom_2nd_v'][1:k]), n['gamma_0'], abs(n['Delta_Sv_hom'])

Ge, GTe, g0, dSv = sweep('elastic')
Gp, GTp, _, _ = sweep('elastoplastic')
G0 = g0 / dSv
GT_el = np.interp(GL_exp, Ge, GTe); GT_ep = np.interp(GL_exp, Gp, GTp)
f_el = (GT_el - G_req) / (GT_el - G0)
f_ep = (GT_ep - G_req) / (GT_ep - G0)

print(f'gamma_0/|dS_V| = {G0:.4e} m K')
print(f'{"tSL[s]":>7}{"P[mm]":>7}{"GL[K/m]":>9}{"SDAS":>7}{"G_req":>11}{"G_el":>11}{"G_EP":>11}{"f_el":>7}{"f_EP":>7}')
for i in range(n_exp):
    print(f'{t_exp[i]:7.1f}{P_exp[i]:7.1f}{GL_exp[i]:9.0f}{s_exp[i]:7.1f}{G_req[i]:11.3e}{GT_el[i]:11.3e}{GT_ep[i]:11.3e}{f_el[i]:7.3f}{f_ep[i]:7.3f}')
m = P_exp <= 90
print('corr(GL, f_EP) (P<=90 mm) =', np.corrcoef(GL_exp[m], f_ep[m])[0, 1])

fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
ax[0].plot(GL_exp[m], G_req[m], 'sk', label=r'$\Gamma$ required by each SDAS point')
ax[0].plot(Gp[(Gp > 900) & (Gp < 11000)], GTp[(Gp > 900) & (Gp < 11000)], '--', color='#eb6834', lw=2, label=r'$\Gamma^{2nd}$ elasto-plastic (loaded)')
ax[0].plot(Ge[(Ge > 900) & (Ge < 11000)], GTe[(Ge > 900) & (Ge < 11000)], '-', color='#2a78d6', lw=2, label=r'$\Gamma^{2nd}$ elastic')
ax[0].axhline(G0, color='#1baf7a', ls='-.', lw=1.5, label=r'$\gamma_0/|\Delta S_V|$ (fully restored)')
ax[0].set_xlabel(r'Local thermal gradient $\nabla T_L$ [K/m]'); ax[0].set_ylabel(r'$\Gamma$ [m K]'); ax[0].legend(fontsize=9)
ax[1].plot(GL_exp[m], f_ep[m], 's', color='#eb6834', label='elasto-plastic surface')
ax[1].plot(GL_exp[m], f_el[m], 'o', color='#2a78d6', mfc='none', label='elastic surface')
ax[1].set_xlabel(r'Local thermal gradient $\nabla T_L$ [K/m]'); ax[1].set_ylabel('Restored fraction of the elastic excess, f'); ax[1].legend(fontsize=9)
for a in ax: a.grid(True, ls='--', color='darkgray', lw=0.5)
fig.suptitle('Al-3Cu-5Nb-0.1Fe: elastic restoration required by each SDAS point (P <= 90 mm)')
fig.tight_layout()
fig.savefig(os.path.join(EP, 'estudo_AlCuNbFe', 'Figure_restauracao_p2p_B_MRS_2026.png'), dpi=150)
