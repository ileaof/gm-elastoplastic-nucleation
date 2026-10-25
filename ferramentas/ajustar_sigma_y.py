"""
Areal yield stress sigma_y of the interface for each alloy, fitted to its SDAS (same procedure as for
Al-3Cu-5Nb-0.1Fe). Only reads the original scripts.

 1. Gamma* : MRB SDAS (prefactor 5.5, Gamma constant) with zero mean log-bias against the experiments
 2. updated reference, H = 0, energy mode of the copies (ENERGY_MODE): sigma_y such that the 2nd-order
    tensor gamma/|dS_V| at the mean gradient equals Gamma* (sweep of gm_elastoplastic)
 3. check that the surface yields below the mean experimental gradient
Writes sigma_y_ligas.txt (read by aplicar_elastoplastico.py).
"""
import os
import re
import sys
import numpy as np
from math import exp, log, pi
from scipy.optimize import brentq

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SRC = os.path.join(ROOT, 'Sim Luane Final')
sys.path.insert(0, SRC)   # midpoint.py and the other helpers of the original scripts
sys.path.insert(0, os.path.join(ROOT, 'modelos_elasto-plastos', 'nucleacao'))
sys.path.insert(0, os.path.dirname(__file__))
from gm_elastoplastic import GMElastoPlasticSurface
from aplicar_elastoplastico import ENERGY_MODE

ALLOYS = {
    # alloy: (SDAS script, nucleation script, mean experimental gradient [K/m])
    'AlCuNbFe': ('scriptsdasmodelAl3Cu5Nb01Fe_zoom_horizontal_B_MRS_Meeting_2026.py',
                 'Nucleation_model_FCC_A1_new_paper_B_MRS_Meeting_2026_dev_full_GLH_AlCuNbFe_signal_beta_2nd_Sim_2026_paper.py', 3294.0),
    'AlSiMgFe': ('scriptsdasmodelAl08Si06Mg02Fe_zoom_horizontal_new_script_Ferreira_2024.py',
                 'Nucleation_model_FCC_A1_new_paper_2024_dev_full_GLH_AlSiMgFe_signal_beta_2nd_Sim_Luane.py', 7902.38),
    'AlCuMgFe': ('scriptsdasmodelAl3Cu05Mg02Fe_zoom_horizontal_new_script_Ferreira_2025.py',
                 'Nucleation_model_FCC_A1_new_paper_MR_2025_dev_full_GLH_AlCuMgFe_signal_beta_2nd_Sim_2026_paper.py', 3708.2),
}


def sdas_namespace(script):
    """Executes the SDAS script up to the first model evaluation (definitions and data only)."""
    src = open(os.path.join(SRC, script), encoding='utf-8').read()
    core = src.split('\nSDAS_RB = calc_RB')[0]
    ns = {}
    cwd = os.getcwd()
    os.chdir(SRC)
    try:
        exec(compile(core, script, 'exec'), ns)
    finally:
        os.chdir(cwd)
    return ns


def constant(src, name):
    m = re.search(r'^%s\s*=\s*([-+0-9.eE]+)' % re.escape(name), src, flags=re.M)
    return float(m.group(1))


lines = ['# alloy  sigma_y[N/m]  Gamma*[m K]  (fitted by ajustar_sigma_y.py)']
for alloy, (sdas, nucl, G_mean) in ALLOYS.items():
    ns = sdas_namespace(sdas)
    P, VL, tSL, xsc, ysc = ns['P'], ns['VL'], ns['tSL'], ns['xsc'], ns['ysc']
    k = (xsc > 0) & (ysc > 0)
    d = ns['diff_length_scale'](P, VL, tSL)

    def bias(GT):
        s, *_ = ns['calc'](GT * np.ones(len(P)), *d)
        return np.mean(np.log(np.interp(xsc[k], tSL, s * 1e6)) - np.log(ysc[k]))

    G_star = brentq(bias, 1e-9, 1e-4, xtol=1e-15)
    s_star, *_ = ns['calc'](G_star * np.ones(len(P)), *d)
    m = np.interp(xsc[k], tSL, s_star * 1e6)
    mape = 100 * np.mean(np.abs(m - ysc[k]) / ysc[k])

    nsrc = open(os.path.join(SRC, nucl), encoding='utf-8').read()
    g0, s0 = constant(nsrc, 'gamma_0'), constant(nsrc, 'sigma_0')
    dSv = constant(nsrc, 'DeltaH') * constant(nsrc, 'rhos') / constant(nsrc, 'TL')
    req = constant(nsrc, 'req')
    gam_star = G_star * dSv
    if gam_star <= g0:
        print(f'{alloy}: Gamma* = {G_star:.4e} needs gamma = {gam_star:.4f} < gamma_0 = {g0}: no sigma_y >= 0')
        continue
    lam, mu_, lam0, mu0_ = (constant(nsrc, n) for n in ('lambda_', 'mu', 'lambda0', 'mu0'))

    def gamma2_at_mean(sy):
        """2nd-order tensor (gamma/|dS_V|) at the mean gradient, sweeping the surface of the copies
        (energy mode ENERGY_MODE, updated reference, H = 0) from r_eq down to grad T = G_mean."""
        os.environ.update(GM_SURFACE_MODEL='elastoplastic', GM_DISSIPATION_MODE=ENERGY_MODE,
                          GM_PLASTIC_REFERENCE='1', GM_H_A='0', GM_SIGMA_YA0=repr(sy))
        surf = GMElastoPlasticSurface(req, g0, s0, lam, mu_, lam0, mu0_)
        Gp = gp = None
        for i in range(1, 20000):
            r = req * (1 - 2.5e-4 * i)
            _, _, g, _ = surf.step(r)
            G = g / (4 * pi * r * r * dSv)
            if G >= G_mean:
                return (gp + (g - gp) * (G_mean - Gp) / (G - Gp)) / dSv if Gp is not None else g / dSv
            Gp, gp = G, g
        raise RuntimeError('mean gradient not reached')

    sy = brentq(lambda x: gamma2_at_mean(x) - G_star, 1e-6, 2 * s0, xtol=1e-7)
    G_y = None
    os.environ['GM_SIGMA_YA0'] = repr(sy)
    surf = GMElastoPlasticSurface(req, g0, s0, lam, mu_, lam0, mu0_)
    for i in range(1, 20000):
        r = req * (1 - 2.5e-4 * i)
        _, _, g, _ = surf.step(r)
        if surf.yielded_at is not None:
            G_y = g / (4 * pi * r * r * dSv)
            break
    print(f'{alloy}: {int(k.sum())} SDAS points, Gamma* = {G_star:.4e} m K (MAPE {mape:.1f}%), '
          f'gamma_0 = {g0}, sigma_0 = {s0}, |dS_V| = {dSv:.4e} -> sigma_y = {sy:.4f} N/m ({ENERGY_MODE}); '
          f'yield at grad T ~ {G_y:.0f} K/m (mean {G_mean:g})')
    lines.append(f'{alloy} {sy:.6g} {G_star:.6e}')

lines.append('H2O 0.0695 nan   # no SDAS data: value of Al-Cu-Nb-Fe kept as a placeholder')
open(os.path.join(os.path.dirname(__file__), 'sigma_y_ligas.txt'), 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
