"""
Elasto-plastic Gurtin-Murdoch surface in non-equilibrium phase nucleation
Solution sequence of the paper "Elasto-Plastic Gurtin-Murdoch Surface in Non-Equilibrium Phase Nucleation:
Stored Permanent Surface Deformation, Elastic Restoration and the Thermal-Field Tensor" (I. L. Ferreira).

Alloy: Al-0.8Si-0.6Mg-0.2Fe (wt.%), FCC_A1 phase, homogeneous nucleus (contact angle theta = pi).
The script solves the elastic and the elasto-plastic surface side by side and writes every output of the paper for
this alloy: yield point, values at the mean experimental gradient and at 10^5 K/m, maximum of the nucleation rate,
the sweep table (CSV) and the figures against the thermal gradient.

ALGORITHM (one radius at a time, r decreasing from the equilibrium radius r_eq)
----------------------------------------------------------------------------------
  0. Input data (Table 2 of the paper)
  for i = 0, 1, 2, ...
  1. Radius                     r_i = r_eq - i * dr,  dr = 0.005 r_eq
  2. Areal strain               E_A = 2 ln[y(r)/y(r_eq)],  y(r) = (3 lambda + 2 mu) r + 2 (lambda_0 + mu_0)    (Sec. 2.3)
                                dE_A/dr = 2 / [r (1 + alpha_r)],  alpha_r = 2 (|sigma_0| + 2 (lambda_0 + mu_0)) / [r (3 lambda + 2 mu)]
  3. Excess surface stress      elastic:        s = C_A E_A,  C_A = |sigma_0|                                  (Sec. 2.3)
                                elasto-plastic: return mapping on |s| <= sigma_y + H alpha                       (Sec. 3.3)
                                  s_trial = C_A (E_A - E_A^p);  phi = |s_trial| - (sigma_y + H alpha)
                                  if phi > 0:  dlam = phi / (C_A + H);  s = s_trial - C_A dlam sign(s_trial)
                                               E_A^p += dlam sign(s_trial);  alpha += dlam
  4. Surface tension            elastic:        tau = sigma_0 + s                                               (Sec. 2.3)
                                elasto-plastic: tau = (sigma_0 + s^e) + |sigma_0| E_A^p + sigma_y C_A/(C_A + H)  (Sec. 3.5)
                                                (elastic component + plastic component + dissipation term)
  5. Surface energy             elastic:        gamma = gamma_0 (r_eq/r)^2 - s                                  (Sec. 2.3)
                                elasto-plastic: gamma = gamma_0 (r_eq/r)^2 + |sigma_0| (|E_A^e| + |E_A^p|) - sigma_y alpha
                                                (permanent deformation stored, dissipated heat removed)          (Sec. 3.5)
     and its radial derivative  d gamma/dr = -2 gamma_0 r_eq^2/r^3 - C_t dE_A/dr - C_A dE_A^p/dr - sigma_y d alpha/dr
                                (consistent tangent C_t = C_A H/(C_A + H) while yielding, C_A otherwise)          (Sec. 3.5)
  6. Undercooling and gradient  Delta T = 2 gamma / (r |Delta S_V|),  grad T = Delta T / (8 pi r)                (Sec. 2.1)
  7. Thermal-field tensors      Gamma_2nd = 4 pi r^2 grad T                                                     (Sec. 2.1)
                                Gamma_1st = -gamma / (Delta S_V + (1/Delta T) d gamma/dr)
  8. Critical free energy       Delta G_C = 4 (pi r^2 gamma + (1/3) pi r^3 Delta S_V Delta T)   (f(theta = pi) = 4)
     and nucleation rate        I = D(T)/a^2 * 4 pi r^2/a^2 * N_V * exp(Delta G_C / Delta G_C,ref)              (Sec. 3.8)
                                (T = T_L - Delta T, Delta G_C,ref = value at the first point of the sweep)
  9. Stop when Delta T > Delta T_max (30 K covers the whole range of solidification gradients, ~10^6 K/m)
 9b. Elastic restoration (Sec. 3.6): only the recoverable elastic energy |s| <= sigma_y is restored, by the fraction
     f = eps_v tau_r / (1 + eps_v tau_r) of the rate of density change; without tau_r the largest restoration, f = 1,
     is shown:  gamma_R = gamma_ep - f |s|,  Gamma_R = Gamma_2nd,ep - f |s| / |Delta S_V|  (> gamma_0/|Delta S_V|)
 9c. Bulk, surface and configurational free energies and entropies (2nd order, Sec. 2.1, Eq. (3)), for i >= 1:
     d(Delta T)/dr from the sweep;  d(Delta S_V)/dr = -3 (Delta S_V Delta T + d gamma/dr)/(2 r Delta T) - (Delta S_V/Delta T) d(Delta T)/dr
     Delta S_V(r_i) = Delta S_V(r_{i-1}) + d(Delta S_V)/dr * dr / Delta T       (accumulated along the sweep)
     1st-order contact angle theta^1st: root of f_het_1st (TERMOA = Delta S_V + d gamma/dr / Delta T,
       TERMOB = gamma, TERMOC = 2 gamma/Delta T), every radius;  2nd-order contact angle theta^2nd: root of
       f_het_2nd;  for both: r_C,het, f(theta), f'(theta) and d theta/dr (central differences along the sweep)
     bulk          Delta G_V = Delta S_V Delta T,          Delta S_V,het = Delta S_V f(theta)/4
     surface       Delta G_S = gamma,                      Delta S_S = (d gamma/dr / Delta T) f(theta)/4
     configurational Delta G_C = gamma f'(theta)/f(theta), Delta S_C = (gamma/Delta T) f'(theta)/f(theta)
     with f(theta) = 2 - 3 cos(theta) + cos^3(theta) and f'(theta) its derivative used in the formulation.
 10. Post-processing: yield point (end of the elastic range), values at the mean experimental gradient and at
     10^5 K/m (linear interpolation in grad T), maximum of the rate, CSV, figures (one file per quantity and the
     panels), comparison with the paper.

Usage:  python nucleation_elastoplastic_AlSiMgFe.py            (needs numpy, scipy and matplotlib)
        options: --show=all (default: every figure on screen with plt.show() after the PNGs), --show=panels
        (only the panel figures and the rate against the undercooling), --show=none or --no-show (PNGs only)
Outputs: in the folder of the script, subfolder outputs_AlSiMgFe/ (individual figures in outputs_AlSiMgFe/figures/)
"""
import csv
import os
from math import cos, exp, log, pi, sin

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import sys
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------------------------ 0. input data
ALLOY = 'Al-0.8Si-0.6Mg-0.2Fe'
T_L = 925.32                 # liquidus temperature [K]
DELTA_H = 335300.0           # latent heat [J/kg]
RHO_S = 2555.72              # solid density [kg/m^3]
DS_V = -DELTA_H * RHO_S / T_L   # volumetric entropy of fusion Delta S_V [J m^-3 K^-1] (negative)
GAMMA_0 = 0.183              # equilibrium surface energy gamma_0 [J/m^2]
SIGMA_0 = 1.09               # |sigma_0|, magnitude of the residual (compressive) surface tension [N/m]
LAMBDA, MU = 5.034e10, 2.594e10      # bulk Lame moduli [Pa]
LAMBDA_0, MU_0 = -5.0, 1.6           # surface Lame moduli [N/m]
R_EQ = 6.3e-6                # equilibrium (reference) radius r_eq = r_0 [m]
DR = 0.005 * R_EQ            # radius step [m]
SIGMA_Y = 0.157175           # areal yield stress of the interface sigma_y^s [N/m]
H = 0.0                      # areal hardening modulus [N/m] (no hardening at the liquidus)
A_JUMP = 0.25e-9             # atomic jump distance a [m]
N_V = 1.1745882323279385e28  # density of modes of the FCC_A1 phase [m^-3] (Debye model of the original script)
GRAD_MEAN = 7902.38          # mean experimental thermal gradient [K/m]
GRAD_MARK = 1.0e5            # gradient at which the values are also reported [K/m]
DT_MAX = 30.0                # end of the sweep [K]
C_A = SIGMA_0                # areal modulus of the excess stress
F_RESTORED = 1.0             # restored fraction of the elastic energy (f = 1: largest elastic restoration)
SHOW = 'all'                 # on screen after the PNGs: 'all', 'panels' or 'none' (command line: --show=panels)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'outputs_AlSiMgFe')

# values of the paper for this alloy (elastic, elasto-plastic), used in the final check
PAPER = {
    'yield gradient [K/m]': (None, 925.0),
    'tau at mean grad [N/m]': (-2.090, -1.980),
    'gamma at mean grad [J/m2]': (1.458, 1.397),
    'Gamma_2nd at mean grad [m K]': (1.575e-6, 1.508e-6),
    'Gamma_1st at mean grad [m K]': (7.637e-7, 7.503e-7),
    'Delta T at mean grad [K]': (0.791, 0.774),
    'r_C at mean grad [um]': (3.982, 3.897),
    'Gamma_2nd at 1e5 K/m [m K]': (4.805e-6, 4.581e-6),
    'gamma at 1e5 K/m [J/m2]': (4.450, 4.242),
    'tau at 1e5 K/m [N/m]': (-3.640, -3.535),
    'I_max [m-3 s-1]': (1.01e50, 8.19e49),
    'grad T at I_max [K/m]': (4465.0, 3632.0),
    'Gamma_2nd restored (f=1), mean grad': (None, None),      # not tabulated in the paper for this alloy
    'Gamma_2nd restored (f=1), 1e5 K/m': (None, None),
    'gamma restored (f=1), mean grad': (None, None),
}


def self_diffusion(T, T_m):
    """Self-diffusion coefficient of liquid Al [m^2/s] (doi:10.1016/j.jnoncrysol.2006.03.049)."""
    D_Tm, p, q = 7.44e-9, 0.10, -0.11
    return D_Tm * (T / T_m) ** 1.5 / ((1 + p - p * T / T_m) * (1 - q + q * T / T_m) ** (2.0 / 3.0))


def f_theta(theta):
    """Shape factor of the spherical cap, f(theta) = 2 - 3 cos(theta) + cos^3(theta) (f(pi) = 4)."""
    return 2 - 3 * cos(theta) + cos(theta) ** 3


def df_theta(theta):
    """Derivative term of f(theta) used in the 2nd-order formulation (as in the original model)."""
    return -3 * (2 - 3 * cos(theta) + cos(theta) ** 3) - (1 - cos(theta)) * (2 - cos(theta) - cos(theta) ** 2)


def r_het(r_hom, theta):
    """Heterogeneous critical radius for the contact angle theta (integral of the original model)."""
    integral, _ = quad(lambda t: sin(t) * (1 + cos(t)) / (2 - cos(t) - cos(t) ** 2), pi, theta)
    return r_hom / exp(-integral)


def f_het_1st(theta, r_hom, DT, gamma, dgamma):
    """1st-order condition for the contact angle (f_het of the original model, TERMOA = Delta S_V + (1/Delta T)
    d gamma/dr, TERMOB = gamma, TERMOC = 2 gamma/Delta T): zero when theta is the 1st-order nucleation angle."""
    a, b, c = DS_V + dgamma / DT, gamma, 2 * gamma / DT
    return -c * f_theta(theta) / 4 / (a * f_theta(theta) / 4 + b * df_theta(theta)) - r_het(r_hom, theta)


def f_het_2nd(theta, r_hom, DT, gamma, dgamma, dDSvDT_dr):
    """2nd-order condition for the contact angle: zero when theta is the 2nd-order nucleation angle."""
    num = (DS_V * DT + dgamma) * f_theta(theta) / 4 + gamma * df_theta(theta)
    den = dDSvDT_dr * f_theta(theta) / 4 + DS_V * DT * df_theta(theta)
    return -1.5 * num / den - r_het(r_hom, theta)


def areal_strain(r):
    """Step 2: areal strain of the spherical Gurtin-Murdoch surface and its radial derivative."""
    y = (3 * LAMBDA + 2 * MU) * r + 2 * (LAMBDA_0 + MU_0)
    y0 = (3 * LAMBDA + 2 * MU) * R_EQ + 2 * (LAMBDA_0 + MU_0)
    alpha_r = 2 * (SIGMA_0 + 2 * (LAMBDA_0 + MU_0)) / (r * (3 * LAMBDA + 2 * MU))
    return 2 * log(y / y0), 2 / (r * (1 + alpha_r))


def solve(plastic):
    """Steps 1-9 for the elastic (plastic=False) or the elasto-plastic (plastic=True) surface.
    Returns a dict of arrays, one value per radius of the sweep."""
    E_p = 0.0          # permanent areal strain E_A^p
    alpha = 0.0        # accumulated areal plastic strain
    rows = []
    i = 0
    while True:
        r = R_EQ - i * DR                                         # 1. radius
        if r <= 0.0:
            break
        E_A, dEA = areal_strain(r)                                # 2. areal strain
        # 3. excess surface stress (return mapping for the elasto-plastic surface)
        s = C_A * (E_A - E_p)
        C_t, dEp, dalpha, yielding = C_A, 0.0, 0.0, False
        if plastic:
            phi = abs(s) - (SIGMA_Y + H * alpha)
            if phi > 0.0:
                dlam = phi / (C_A + H)
                sgn = 1.0 if s > 0.0 else -1.0
                s -= C_A * dlam * sgn
                E_p += dlam * sgn
                alpha += dlam
                C_t = C_A * H / (C_A + H)                         # consistent tangent
                dEp = C_A / (C_A + H) * dEA                       # dE_A^p/dr
                dalpha = sgn * C_A / (C_A + H) * dEA              # d alpha/dr
                yielding = True
        # 4. surface tension: elastic component, plastic component and dissipation term
        tau_e = -SIGMA_0 + s                                      # sigma_0 + s^e (sigma_0 < 0)
        tau_p = C_A * E_p                                         # |sigma_0| E_A^p
        tau_d = SIGMA_Y * C_A / (C_A + H) if alpha > 0.0 else 0.0
        tau = tau_e + tau_p + tau_d
        # 5. surface energy and its radial derivative
        gamma = GAMMA_0 * (R_EQ / r) ** 2 - s - C_A * E_p - SIGMA_Y * alpha
        dgamma = -2 * GAMMA_0 * R_EQ ** 2 / r ** 3 - C_t * dEA - C_A * dEp - SIGMA_Y * dalpha
        # 6. undercooling and thermal gradient (closed form of the critical condition)
        DT = 2 * gamma / (r * abs(DS_V))
        if DT > DT_MAX:                                           # 9. end of the sweep
            break
        grad = DT / (8 * pi * r)
        # 7. thermal-field tensors
        G2 = 4 * pi * r ** 2 * grad
        G1 = -gamma / (DS_V + dgamma / DT)
        # 8. critical free energy (homogeneous, f(pi) = 4) and the nucleation-rate prefactor
        DGc = 4 * (pi * r ** 2 * gamma + pi * r ** 3 * DS_V * DT / 3.0)
        rows.append(dict(r=r, E_A=E_A, E_p=E_p, alpha=alpha, s=s, tau=tau, tau_e=tau_e, tau_p=tau_p, gamma=gamma,
                         dgamma=dgamma, DT=DT, grad=grad, G2=G2, G1=G1, DGc=DGc, yielding=yielding))
        i += 1
    out = {k: np.array([row[k] for row in rows]) for k in rows[0]}
    # 9b. elastic restoration, largest (f = 1): only the elastic part |s| of the surface energy is restored
    out['gamma_R'] = out['gamma'] - F_RESTORED * np.abs(out['s'])
    out['G2_R'] = out['G2'] - F_RESTORED * np.abs(out['s']) / abs(DS_V)
    # 8. nucleation rate, referred to the critical free energy of the first point after r_eq (index 1)
    T = T_L - out['DT']
    out['I'] = (self_diffusion(T, T_L) / A_JUMP ** 2) * (4 * pi * out['r'] ** 2 / A_JUMP ** 2) * N_V \
        * np.exp(out['DGc'] / out['DGc'][1])
    # 9c. bulk, surface and configurational free energies and entropies
    second_order_energies(out)
    return out


def second_order_energies(out):
    """Step 9c: 2nd-order contact angle, bulk, surface and configurational free energies and entropies."""
    n = len(out['r'])
    for k in ('theta', 'r_het', 'f2', 'df2', 'theta1', 'r_het1', 'f1', 'df1',
              'DSv', 'G_V', 'G_S', 'G_C', 'S_V', 'S_S', 'S_C'):
        out[k] = np.full(n, np.nan)
    out['DSv'][0] = DS_V
    for i in range(n):                     # 1st-order contact angle: root of f_het_1st (every radius)
        r, DT, gamma, dgamma = out['r'][i], out['DT'][i], out['gamma'][i], out['dgamma'][i]
        try:
            theta1 = brentq(f_het_1st, 0.01, pi, args=(r, DT, gamma, dgamma), xtol=2e-14)
        except (ValueError, RuntimeError):
            continue
        out['theta1'][i], out['r_het1'][i] = theta1, r_het(r, theta1)
        out['f1'][i], out['df1'][i] = f_theta(theta1), df_theta(theta1)
    for i in range(1, n):
        r, DT, gamma, dgamma = out['r'][i], out['DT'][i], out['gamma'][i], out['dgamma'][i]
        dDT_dr = (DT - out['DT'][i - 1]) / (r - out['r'][i - 1])
        dDSv_dr = -3 * (DS_V * DT + dgamma) / (2 * r * DT) - DS_V / DT * dDT_dr
        out['DSv'][i] = out['DSv'][i - 1] + dDSv_dr * DR / DT
        dDSvDT_dr = -1.5 * (DS_V * DT + dgamma) / r
        try:
            theta = brentq(f_het_2nd, 0.01, pi, args=(r, DT, gamma, dgamma, dDSvDT_dr), xtol=2e-14)
        except ValueError:                 # no 2nd-order angle: the heterogeneous quantities end here
            break
        f, df = f_theta(theta), df_theta(theta)
        out['theta'][i], out['r_het'][i] = theta, r_het(r, theta)
        out['f2'][i], out['df2'][i] = f, df
        out['G_V'][i] = out['DSv'][i] * DT                 # bulk free energy        [J/m^3]
        out['G_S'][i] = gamma                              # surface free energy     [J/m^2]
        out['G_C'][i] = gamma / f * df                     # configurational energy  [J/m^2]
        out['S_V'][i] = out['DSv'][i] * f / 4              # bulk entropy            [J m^-3 K^-1]
        out['S_S'][i] = dgamma / DT * f / 4                # surface entropy         [J m^-3 K^-1]
        out['S_C'][i] = gamma / DT / f * df                # configurational entropy [J m^-2 K^-1]
    # angles in degrees and their radial derivatives d theta/dr (central differences along the sweep)
    for o in ('1', '2'):
        th = out['theta1' if o == '1' else 'theta']
        out[f'theta{o}_deg'] = np.degrees(th)
        out[f'dtheta{o}_dr'] = np.gradient(th, out['r'])


def at(sol, key, grad):
    """Value of a quantity at a given thermal gradient (linear interpolation along the sweep)."""
    k = slice(1, None)
    return float(np.interp(grad, sol['grad'][k], sol[key][k]))


def main():
    os.makedirs(OUT, exist_ok=True)
    el, ep = solve(False), solve(True)

    # ------------------------------------------------------------------------------------ 10. post-processing
    iy = int(np.argmax(ep['yielding']))
    res = {
        'yield gradient [K/m]': (None, ep['grad'][iy]),
        'tau at mean grad [N/m]': (at(el, 'tau', GRAD_MEAN), at(ep, 'tau', GRAD_MEAN)),
        'gamma at mean grad [J/m2]': (at(el, 'gamma', GRAD_MEAN), at(ep, 'gamma', GRAD_MEAN)),
        'Gamma_2nd at mean grad [m K]': (at(el, 'G2', GRAD_MEAN), at(ep, 'G2', GRAD_MEAN)),
        'Gamma_1st at mean grad [m K]': (at(el, 'G1', GRAD_MEAN), at(ep, 'G1', GRAD_MEAN)),
        'Delta T at mean grad [K]': (at(el, 'DT', GRAD_MEAN), at(ep, 'DT', GRAD_MEAN)),
        'r_C at mean grad [um]': (1e6 * at(el, 'r', GRAD_MEAN), 1e6 * at(ep, 'r', GRAD_MEAN)),
        'Gamma_2nd at 1e5 K/m [m K]': (at(el, 'G2', GRAD_MARK), at(ep, 'G2', GRAD_MARK)),
        'gamma at 1e5 K/m [J/m2]': (at(el, 'gamma', GRAD_MARK), at(ep, 'gamma', GRAD_MARK)),
        'tau at 1e5 K/m [N/m]': (at(el, 'tau', GRAD_MARK), at(ep, 'tau', GRAD_MARK)),
        'Gamma_2nd restored (f=1), mean grad': (None, at(ep, 'G2_R', GRAD_MEAN)),
        'Gamma_2nd restored (f=1), 1e5 K/m': (None, at(ep, 'G2_R', GRAD_MARK)),
        'gamma restored (f=1), mean grad': (None, at(ep, 'gamma_R', GRAD_MEAN)),
        'I_max [m-3 s-1]': (el['I'][1:].max(), ep['I'][1:].max()),
        'grad T at I_max [K/m]': (el['grad'][1:][np.argmax(el['I'][1:])], ep['grad'][1:][np.argmax(ep['I'][1:])]),
    }
    print(f'{ALLOY}, FCC_A1, homogeneous nucleus;  sigma_y = {SIGMA_Y} N/m, H = {H} N/m')
    print(f'yield (end of the elastic range): r/r_eq = {ep["r"][iy] / R_EQ:.3f}, Delta T = {ep["DT"][iy]:.3f} K, '
          f'grad T = {ep["grad"][iy]:.0f} K/m')
    print(f'{"quantity":36s}{"elastic":>10s}{"elasto-plastic":>16s}{"paper (el / ep)":>26s}  check')
    for k, (ve, vp) in res.items():
        pe, pp = PAPER[k]
        ok = all(p is None or abs(v / p - 1) < 0.01 for v, p in ((ve, pe), (vp, pp)))
        fmt = lambda v: '' if v is None else f'{v:.4g}'
        check = '' if (pe is None and pp is None) else ('ok' if ok else 'DIFFERS')
        print(f'{k:36s}{fmt(ve):>10s}{fmt(vp):>16s}{fmt(pe) + " / " + fmt(pp):>26s}  {check}')
    print(f'equilibrium: gamma_0 = {GAMMA_0} J/m2, gamma_0/|Delta S_V| = {GAMMA_0 / abs(DS_V):.4e} m K '
          f'(the restored tensor stays above it)')

    # bulk, surface and configurational free energies and entropies at the mean experimental gradient
    print('contact angles (1st and 2nd order), free energies and entropies at the mean gradient '
          '(elastic | elasto-plastic):')
    for key, name in (('theta1_deg', '1st-order contact angle theta [deg]'), ('f1', '1st-order f(theta)'),
                      ('df1', "1st-order f'(theta)"), ('dtheta1_dr', '1st-order d theta/dr [rad/m]'),
                      ('r_het1', '1st-order r_C,het [m]'),
                      ('theta2_deg', '2nd-order contact angle theta [deg]'), ('f2', '2nd-order f(theta)'),
                      ('df2', "2nd-order f'(theta)"), ('dtheta2_dr', '2nd-order d theta/dr [rad/m]'),
                      ('r_het', '2nd-order r_C,het [m]'), ('G_V', 'bulk free energy Delta G_V [J/m3]'),
                      ('G_S', 'surface free energy Delta G_S [J/m2]'), ('G_C', 'configurational Delta G_C [J/m2]'),
                      ('S_V', 'bulk entropy Delta S_V [J/(m3 K)]'), ('S_S', 'surface entropy Delta S_S [J/(m3 K)]'),
                      ('S_C', 'configurational entropy Delta S_C [J/(m2 K)]')):
        print(f'   {name:46s}{at(el, key, GRAD_MEAN):>14.5g}{at(ep, key, GRAD_MEAN):>14.5g}')

    # sweep table
    with open(os.path.join(OUT, 'sweep_AlSiMgFe.csv'), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter=';')
        keys = ['r', 'grad', 'DT', 'tau', 'tau_e', 'tau_p', 'gamma', 'gamma_R', 'dgamma', 'G1', 'G2', 'G2_R', 'DGc', 'I',
                'theta1', 'f1', 'df1', 'dtheta1_dr', 'r_het1', 'theta', 'f2', 'df2', 'dtheta2_dr', 'r_het', 'G_V', 'G_S', 'G_C', 'S_V', 'S_S', 'S_C', 'E_A', 'E_p', 'alpha']
        w.writerow(['surface'] + keys)
        for name, sol in (('elastic', el), ('elasto-plastic', ep)):
            for j in range(len(sol['r'])):
                w.writerow([name] + [f'{sol[k][j]:.6e}' for k in keys])

    # ------------------------------------------------------------------------------------------------ figures
    G_y = ep['grad'][iy]
    show = next((a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith('--show=')), SHOW)
    if '--no-show' in sys.argv[1:]:
        show = 'none'

    def draw(ax, key, label, logy):
        """One quantity against the thermal gradient: elastic (solid), elasto-plastic (dashed), elastic range
        (shaded), yield (dash-dot), mean experimental gradient (dotted) and 10^5 K/m (red)."""
        k = np.isfinite(el[key]) & (np.arange(len(el[key])) > 0)
        kp = np.isfinite(ep[key]) & (np.arange(len(ep[key])) > 0)
        ax.axvspan(el['grad'][1], G_y, color='#2a78d6', alpha=0.08, lw=0)
        ax.plot(el['grad'][k], el[key][k], color='#2a78d6', lw=2, label='Elastic GM surface')
        ax.plot(ep['grad'][kp], ep[key][kp], color='#eb6834', lw=2, ls='--', label='Elasto-plastic GM surface')
        if key == 'tau':
            ax.plot(ep['grad'][1:], ep['tau_e'][1:], color='k', lw=1, ls=':', label=r'Elastic component $\sigma_0+s^e$')
            ax.plot(ep['grad'][1:], ep['tau_p'][1:], color='#8a3ffc', lw=1, ls='-.', label=r'Plastic component $|\sigma_0|E_A^p$')
        if key in ('gamma', 'G2'):     # 9b. elastic restoration (f = 1) and equilibrium value
            ax.plot(ep['grad'][1:], ep[key + '_R'][1:], color='k', lw=1.3, ls=':',
                    label='Elasto-plastic, elastic part restored ($f$ = 1)')
            ax.axhline(GAMMA_0 if key == 'gamma' else GAMMA_0 / abs(DS_V), color='0.5', lw=1, ls='-',
                       label='Equilibrium ' + (r'$\gamma_0$' if key == 'gamma' else r'$\gamma_0/|\Delta S_V|$'))
        ax.axvline(G_y, color='k', lw=1, ls='-.')
        ax.axvline(GRAD_MEAN, color='0.4', lw=1, ls=':')
        ax.axvline(GRAD_MARK, color='#b02a2a', lw=1, ls='--')
        ax.set_xscale('log')
        if logy:
            ax.set_yscale('log')
        else:
            ax.ticklabel_format(axis='y', useOffset=False)
        ax.set_xlabel(r'Thermal gradient $\nabla T$ [K m$^{-1}$]')
        ax.set_ylabel(label)
        ax.grid(True, which='both', ls='--', color='0.85', lw=0.5)
        ax.legend(fontsize=7.5, loc='best')

    nucleation = [('tau', r'Surface tension $\tau_{SL}$ [N m$^{-1}$]', False),
                  ('gamma', r'Surface energy $\gamma_{SL}$ [J m$^{-2}$]', True),
                  ('dgamma', r'$\partial\gamma_{SL}/\partial r$ [J m$^{-3}$]', False),
                  ('G2', r'2nd-order thermal-field tensor $\Gamma^{2nd}$ [m K]', True),
                  ('G1', r'1st-order thermal-field tensor $\Gamma^{1st}$ [m K]', True),
                  ('DT', r'Undercooling $\Delta T$ [K]', True),
                  ('r', r'Critical radius $r_C$ [m]', False),
                  ('DGc', r'Critical free energy $\Delta G_C$ [J]', False),
                  ('I', r'Nucleation rate $I$ [m$^{-3}$ s$^{-1}$]', False)]
    energies = [('G_V', r'Bulk free energy $\Delta G_V=\Delta S_V\Delta T$ [J m$^{-3}$]', False),
                ('G_S', r'Surface free energy $\Delta G_S$ [J m$^{-2}$]', True),
                ('G_C', r'Configurational free energy $\Delta G_C$ [J m$^{-2}$]', False),
                ('S_V', r'Bulk entropy $\Delta S_V$ [J m$^{-3}$ K$^{-1}$]', False),
                ('S_S', r'Surface entropy $\Delta S_S$ [J m$^{-3}$ K$^{-1}$]', False),
                ('S_C', r'Configurational entropy $\Delta S_C$ [J m$^{-2}$ K$^{-1}$]', False)]
    angles = [('theta1_deg', r'1st-order contact angle $\theta^{1st}$ [deg]', False),
              ('f1', r'1st-order $f(\theta^{1st})=2-3\cos\theta+\cos^3\theta$', False),
              ('df1', r"1st-order $f'(\theta^{1st})$", False),
              ('dtheta1_dr', r'1st-order $\partial\theta^{1st}/\partial r$ [rad m$^{-1}$]', False),
              ('r_het1', r'1st-order heterogeneous radius $r_{C,het}^{1st}$ [m]', False),
              ('theta2_deg', r'2nd-order contact angle $\theta^{2nd}$ [deg]', False),
              ('f2', r'2nd-order $f(\theta^{2nd})=2-3\cos\theta+\cos^3\theta$', False),
              ('df2', r"2nd-order $f'(\theta^{2nd})$", False),
              ('dtheta2_dr', r'2nd-order $\partial\theta^{2nd}/\partial r$ [rad m$^{-1}$]', False),
              ('r_het', r'2nd-order heterogeneous radius $r_{C,het}^{2nd}$ [m]', False)]

    # one file per quantity
    fig_dir = os.path.join(OUT, 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    for key, label, logy in nucleation + energies + angles:
        fig, ax = plt.subplots(figsize=(7, 5))
        draw(ax, key, label, logy)
        ax.set_title(ALLOY, fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, f'{key}_vs_gradient_AlSiMgFe.png'), dpi=150)
        if show != 'all':
            plt.close(fig)

    # the same figures as panels
    for name, group, shape in (('curves_vs_gradient', nucleation, (3, 3)),
                               ('energies_entropies_vs_gradient', energies, (2, 3)),
                               ('contact_angles_vs_gradient', angles, (2, 5))):
        fig, axs = plt.subplots(*shape, figsize=(5 * shape[1], 4.2 * shape[0]))
        for ax, (key, label, logy) in zip(axs.flat, group):
            draw(ax, key, label, logy)
        fig.suptitle(f'{ALLOY}: elastic and elasto-plastic Gurtin-Murdoch surface (shaded: elastic range; '
                     f'dash-dot: yield; dotted: mean gradient; red: 10$^5$ K/m)', fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        fig.savefig(os.path.join(OUT, f'{name}_AlSiMgFe.png'), dpi=150)
        if show == 'none':
            plt.close(fig)

    # nucleation rate against the undercooling (as in the original scripts)
    fig, ax = plt.subplots(figsize=(7, 4.8))
    ax.plot(el['DT'][1:], el['I'][1:], color='#2a78d6', lw=2, label='Elastic GM surface')
    ax.plot(ep['DT'][1:], ep['I'][1:], color='#eb6834', lw=2, ls='--', label='Elasto-plastic GM surface')
    ax.axvline(ep['DT'][iy], color='k', lw=1, ls='-.', label='Yield (end of the elastic range)')
    ax.set_xlabel(r'Undercooling $\Delta T$ [K]')
    ax.set_ylabel(r'Nucleation rate $I$ [m$^{-3}$ s$^{-1}$]')
    ax.set_title(f'{ALLOY}: homogeneous nucleation rate')
    ax.grid(True, ls='--', color='0.85', lw=0.5)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'rate_vs_undercooling_AlSiMgFe.png'), dpi=150)
    fig.savefig(os.path.join(fig_dir, 'I_vs_undercooling_AlSiMgFe.png'), dpi=150)
    if show == 'none':
        plt.close(fig)
    print(f'outputs written to {OUT}')
    # on screen, after every PNG has been written
    if show != 'none':
        plt.show()


if __name__ == '__main__':
    main()
