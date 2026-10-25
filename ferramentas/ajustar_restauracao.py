"""
Elastic restoration of the surface driven by the rate of density change (Al-3Cu-5Nb-0.1Fe, the alloy with
the local thermal gradient G_L(P) measured along the casting). Only reads the original scripts.

  rate of density change   eps_v = (1/rho) d rho/dt = beta_s * f_s_dot = beta_s * T_dot / DT_SL
                           beta_s = (rho_S - rho_L) / rho_S,   T_dot = G_L(P) * V_L(P)
  solidification interval  DT_SL = 1/(90 - 5) * integral_5^90 T_dot(P) * t_SL(P) dP   (mean integral, measured range)
  restored fraction        f = eps_v tau_r / (1 + eps_v tau_r)
  2nd-order tensor         Gamma(P) = Gamma_EP(G_L) - f * |s(G_L)| / |DS_V|
                           only the elastic (recoverable) energy |s| <= sigma_y is restored; the stored energy of
                           the permanent deformation stays (author, 2026-10-22): Gamma > gamma_0/|DS_V| always
                           Gamma_EP: 2nd-order heterogeneous tensor of the elasto-plastic copy (loaded surface)
  tau_r: environment variable RESTORATION_TAU [s]; without it, f = 1 (largest elastic restoration).
G_L(P) is the quadratic fit 5-90 mm of Thermal_data_Ivaldo_B_MRS_Meeting_2026_Al3Cu5Nb_GL_data_indev.py,
held constant outside. Writes restauracao_AlCuNbFe.txt and gamma_restaurado_AlCuNbFe.csv (read by aplicar_sdas.py).
"""
import os
import sys
import runpy
import tempfile
import numpy as np
from scipy import integrate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SRC = os.path.join(ROOT, 'Sim Luane Final')
NUC = os.path.join(ROOT, 'modelos_elasto-plastos', 'nucleacao')
sys.path.insert(0, SRC)
sys.path.insert(0, NUC)

SDAS_SCRIPT = 'scriptsdasmodelAl3Cu5Nb01Fe_zoom_horizontal_B_MRS_Meeting_2026.py'
NUCL_EP = 'Nucleation_model_FCC_A1_new_paper_B_MRS_Meeting_2026_dev_full_GLH_AlCuNbFe_signal_beta_2nd_Sim_2026_paper_EP.py'
P1, P2 = 5.0, 90.0   # measured range of G_L [mm]
N_EXP = 16


def grad_exp_fit(P):
    """Quadratic fit of G_L(P) [K/m] (P in mm), Thermal_data_Ivaldo_B_MRS_Meeting_2026_Al3Cu5Nb_GL_data_indev.py."""
    P = np.clip(P, P1, P2)
    return 10405.54 - 200.22 * P + 1.04692 * P ** 2


# SDAS model of the original script (definitions and data only)
src = open(os.path.join(SRC, SDAS_SCRIPT), encoding='utf-8').read()
ns = {}
cwd = os.getcwd()
os.chdir(SRC)
try:
    exec(compile(src.split('\nSDAS_RB = calc_RB')[0], SDAS_SCRIPT, 'exec'), ns)
finally:
    os.chdir(cwd)
P, VL, tSL, xsc, ysc = ns['P'], ns['VL'], ns['tSL'], ns['xsc'][:N_EXP], ns['ysc'][:N_EXP]
d = ns['diff_length_scale'](P, VL, tSL)
beta_s = (ns['rhos'] - ns['rhol']) / ns['rhos']

# elasto-plastic nucleation sweep (loaded surface, default mode of the copies)
os.environ['GM_SURFACE_MODEL'] = 'elastoplastic'
plt.show = lambda *a, **k: plt.close('all')
with tempfile.TemporaryDirectory() as wd:
    os.chdir(wd)
    try:
        nn = runpy.run_path(os.path.join(NUC, NUCL_EP), run_name='__main__')
    finally:
        os.chdir(cwd)
        plt.close('all')
k = nn['npoints']
G_sw, GT_sw = np.asarray(nn['dTdr_v'][1:k]), np.asarray(nn['GT_het_2nd_v'][1:k])
s_sw = np.abs(np.asarray(nn['surface_stress_v'][1:k]))   # |s|: recoverable (elastic) excess surface energy
DSv = abs(nn['Delta_Sv_hom'])
GT_0 = nn['gamma_0'] / abs(nn['Delta_Sv_hom'])

# rate of density change
GL = grad_exp_fit(P)
Tdot = GL * VL * 1e-3                                   # V_L [mm/s]
m = (P >= P1) & (P <= P2)
DT_SL = integrate.simpson(Tdot[m] * tSL[m], x=P[m]) / (P2 - P1)
eps_v = beta_s * Tdot / DT_SL
GT_EP = np.interp(GL, G_sw, GT_sw)
GT_el = np.interp(GL, G_sw, s_sw) / DSv               # restorable part of the tensor


def gamma_restored(tau):
    f = np.ones_like(eps_v) if np.isinf(tau) else eps_v * tau / (1 + eps_v * tau)
    return f, GT_EP - f * GT_el


def sdas(GT):
    s, *_ = ns['calc'](GT, *d)
    return s * 1e6


P_exp = (xsc / 4.93116) ** (1 / 0.87606)               # tSL = 4.93116 P^0.87606 of the SDAS script
use = P_exp <= P2


tau_r = float(os.environ.get('RESTORATION_TAU', 'inf'))
f, GT_rest = gamma_restored(tau_r)
S = np.interp(xsc, tSL, sdas(GT_rest))
S_ep = np.interp(xsc, tSL, sdas(GT_EP))
mape = lambda s: 100 * np.mean(np.abs(s[use] - ysc[use]) / ysc[use])

print(f'beta_s = {beta_s:.4f}, DT_SL (mean integral 5-90 mm) = {DT_SL:.2f} K, tau_r = {tau_r:.4g} s, '
      f'gamma_0/|DS_V| = {GT_0:.4e} m K')
print(f'MAPE (P <= 90 mm): loaded elasto-plastic {mape(S_ep):.1f} %, with elastic restoration {mape(S):.1f} %; '
      f'largest restoration of Gamma = {GT_el.max():.3e} m K, Gamma/Gamma_0 >= {(GT_rest/GT_0).min():.2f}')
print(f'{"P[mm]":>6}{"G_L":>8}{"T_dot":>8}{"eps_v":>11}{"f":>7}{"Gamma_EP":>11}{"Gamma":>11}{"SDAS":>8}{"exp":>7}')
for i in range(N_EXP):
    j = int(np.argmin(np.abs(tSL - xsc[i])))
    print(f'{P_exp[i]:6.1f}{GL[j]:8.0f}{Tdot[j]:8.3f}{eps_v[j]:11.3e}{f[j]:7.3f}{GT_EP[j]:11.3e}{GT_rest[j]:11.3e}{S[i]:8.1f}{ysc[i]:7.1f}')

with open(os.path.join(HERE, 'restauracao_AlCuNbFe.txt'), 'w') as fh:
    fh.write('# elastic restoration driven by the rate of density change (ajustar_restauracao.py)\n')
    fh.write(f'beta_s {beta_s:.6g}\nDT_SL {DT_SL:.6g}\ntau_r {tau_r:.6g}\nGamma_0 {GT_0:.6e}\n')
    # restorable part of the tensor at the mean experimental gradient (constant-Gamma curve of the SDAS copy)
    fh.write(f'dGamma_el_mean {np.interp(3294.0, G_sw, s_sw) / DSv:.6e}\n')
    fh.write(f'MAPE_loaded {mape(S_ep):.3f}\nMAPE_restored {mape(S):.3f}\n')
with open(os.path.join(HERE, 'gamma_restaurado_AlCuNbFe.csv'), 'w') as fh:
    fh.write('P_mm;GL;Tdot;eps_v;f;GT_EP;GT_restored\n')
    for i in range(len(P)):
        fh.write(f'{P[i]:g};{GL[i]:.6g};{Tdot[i]:.6g};{eps_v[i]:.6e};{f[i]:.6f};{GT_EP[i]:.6e};{GT_rest[i]:.6e}\n')
