"""Butler equation of the solid-liquid interface of the FCC_A1 nucleus of Al-0.8Si-0.6Mg-0.2Fe (wt.%).

Standard Butler (no equilibrium shift):
    gamma = gamma_i + RT/A_i ln(X_i^S/X_i^B) + (1/A_i)[G_i^{E,S}(T,X^S) - G_i^{E,B}(T,X^B)],   i = Al, Si, Mg, Fe
at the liquidus T_L.  Butler with equilibrium shift (Ferreira, Materials 2026, 19, 4676, Eq. 34), fed by the
nucleation data of the elastic and elasto-plastic Gurtin-Murdoch surfaces (outputs_AlSiMgFe/sweep_AlSiMgFe.csv):
    T = T_L - Delta T(grad T),    G_i^{E,C}/A_i = X_i^S Delta G_C(grad T),   Delta G_C = gamma f'(theta)/f(theta)
so the configurational (substrate) term enters each component in proportion to its occupation of the surface.
Second series (start of Butler = start of nucleation): the configurational term counted from equilibrium,
    Delta G_C^rel = (f'/f)(gamma - gamma_0) = Delta G_C (1 - gamma_0/gamma),
which vanishes at the equilibrium surface energy gamma_0, so Butler with shift starts from standard Butler.
Unknowns X_Si^S, X_Mg^S, X_Fe^S and gamma (X_Al^S = 1 - sum): four Butler equations.

Liquid excess Gibbs energy: Redlich-Kister + Muggianu with the COST507 liquid parameters (COST507-modified.tdb,
gist bocklund/c4714ddbc0500c78e6fe255a763e7550); surface excess G^{E,S} = beta G^{E,B}(T, X^S), beta = 0.83 (Tanaka).
Pure-component interfacial energies by Turnbull's rule gamma_i = alpha_T dH_m,i/(N_A^{1/3} V_i^{2/3}) and molar
interface areas A_i = 1.091 N_A^{1/3} V_i^{2/3}.
"""
import csv
import os
from math import log

import numpy as np
from scipy.optimize import brentq, fsolve

HERE = os.path.dirname(os.path.abspath(__file__))
SWEEP = os.path.join(HERE, 'outputs_AlSiMgFe', 'sweep_AlSiMgFe.csv')
OUT = os.path.join(HERE, 'outputs_AlSiMgFe')

R, N_A = 8.314462618, 6.02214076e23
T_L = 925.32                                  # liquidus [K] (nucleation model)
GAMMA_0 = 0.183                               # equilibrium surface energy of the nucleus [J/m^2] (nucleation model)
EL = ('AL', 'SI', 'MG', 'FE')
WT = {'SI': 0.8, 'MG': 0.6, 'FE': 0.2}        # wt.%
M = {'AL': 26.9815, 'SI': 28.0855, 'MG': 24.305, 'FE': 55.845}       # g/mol
V = {'AL': 10.00e-6, 'SI': 12.06e-6, 'MG': 14.00e-6, 'FE': 7.09e-6}  # molar volumes of the solids [m^3/mol]
DH_M = {'AL': 10710.0, 'SI': 50210.0, 'MG': 8480.0, 'FE': 13810.0}   # enthalpies of fusion [J/mol]
ALPHA_T = {'AL': 0.45, 'SI': 0.32, 'MG': 0.45, 'FE': 0.45}           # Turnbull: metals 0.45, semimetals 0.32
F_AREA = 1.091                                # geometric factor of the close-packed interface
BETA = 0.83                                   # surface/bulk ratio of the excess Gibbs energy

# COST507 liquid interaction parameters L_k(T) = a + b T [J/mol] (binary i-j: (x_i - x_j)^k)
L_BIN = {
    ('AL', 'FE'): [(-91976.5, 22.1314), (-5672.58, 4.8728), (121.9, 0.0)],
    ('AL', 'MG'): [(-12000.0, 8.566), (1894.0, -3.000), (2000.0, 0.0)],
    ('AL', 'SI'): [(-11340.10, -1.23394), (-3530.93, 1.35993), (2265.39, 0.0)],
    ('FE', 'MG'): [(61343.0, 1.5), (-2700.0, 0.0)],
    ('FE', 'SI'): [(-164434.6, 41.9773), (0.0, -21.523), (-18821.542, 22.07), (9695.8, 0.0)],
    ('MG', 'SI'): [(-83864.26, 32.44438), (18027.41, -19.61202), (2486.67, -0.31084),
                   (18541.17, -2.31766), (-12338.84, 1.54236)],
}
L_TER = {('AL', 'MG', 'SI'): [(4125.86, -0.51573), (-47961.64, 5.99520), (25813.80, -3.22672)]}
IDX = {e: k for k, e in enumerate(EL)}

A = {e: F_AREA * N_A ** (1 / 3) * V[e] ** (2 / 3) for e in EL}                          # [m^2/mol]
GAMMA_I = {e: ALPHA_T[e] * DH_M[e] / (N_A ** (1 / 3) * V[e] ** (2 / 3)) for e in EL}     # [J/m^2]


def bulk_fractions():
    n = {e: WT[e] / M[e] for e in WT}
    n['AL'] = (100 - sum(WT.values())) / M['AL']
    tot = sum(n.values())
    return np.array([n[e] / tot for e in EL])


def g_excess(x, T):
    """Molar excess Gibbs energy of the liquid [J/mol] (Redlich-Kister, Muggianu)."""
    g = 0.0
    for (i, j), Ls in L_BIN.items():
        xi, xj = x[IDX[i]], x[IDX[j]]
        g += xi * xj * sum((a + b * T) * (xi - xj) ** k for k, (a, b) in enumerate(Ls))
    for (i, j, k), Ls in L_TER.items():
        xi, xj, xk = x[IDX[i]], x[IDX[j]], x[IDX[k]]
        d = (1 - xi - xj - xk) / 3
        g += xi * xj * xk * sum((a + b * T) * v for (a, b), v in zip(Ls, (xi + d, xj + d, xk + d)))
    return g


def partial_excess(x, T, h=1e-7):
    """Partial molar excess Gibbs energies G_i^E = G^E + dG^E/dx_i - sum_j x_j dG^E/dx_j [J/mol]."""
    grad = np.empty(len(x))
    for k in range(len(x)):
        xp, xm = x.copy(), x.copy()
        xp[k] += h
        xm[k] -= h
        grad[k] = (g_excess(xp, T) - g_excess(xm, T)) / (2 * h)
    return g_excess(x, T) + grad - x @ grad


def butler(T, xb, dGc=0.0, guess=None, phi=None):
    """Solve the Butler system; dGc is the configurational free energy Delta G_C [J/m^2] (0: standard Butler).
    phi: optional element-specific term phi_i [J/m^2], entered as phi_i - sum_j X_j^S phi_j (zero mean)."""
    gb = partial_excess(xb, T)
    a = np.array([A[e] for e in EL])
    g0 = np.array([GAMMA_I[e] for e in EL])
    phi = np.zeros(4) if phi is None else np.asarray(phi, float)
    pbar = float(xb @ phi)

    def xs_given(gam, gs):
        # each Butler equation is monotonic in its own X_i^S: gamma_i + RT/A_i ln(x/x_B) - x dGc + ... = gamma
        xs = np.empty(4)
        for k in range(4):
            f = lambda lx: g0[k] + R * T / a[k] * (lx - log(xb[k])) + (gs[k] - gb[k]) / a[k] - np.exp(lx) * dGc \
                + phi[k] - pbar - gam
            if f(0.0) <= 0:
                xs[k] = 1.0
            elif f(-200.0) >= 0:
                xs[k] = 0.0
            else:
                xs[k] = np.exp(brentq(f, -200.0, 0.0, xtol=1e-15))
        return xs

    # initial guess: damped fixed point on the surface excess energies, each step solved exactly
    gs = BETA * partial_excess(xb if guess is None else guess, T)
    for _ in range(40):
        gam = brentq(lambda g: xs_given(g, gs).sum() - 1.0, -50.0, 50.0 + abs(dGc), xtol=1e-15)
        xs = xs_given(gam, gs)
        gs = 0.5 * gs + 0.5 * BETA * partial_excess(xs / xs.sum(), T)
        pbar = 0.5 * pbar + 0.5 * float(xs @ phi) / xs.sum()

    # Newton polish of the coupled system (unknowns ln X_Si^S, ln X_Mg^S, ln X_Fe^S, gamma)
    def butler_res(u):
        xs = np.concatenate(([1 - np.exp(u[:3]).sum()], np.exp(u[:3])))
        if xs[0] <= 0:
            return np.full(4, 1e3)
        return g0 + R * T / a * np.log(xs / xb) + (BETA * partial_excess(xs, T) - gb) / a - xs * dGc + phi - xs @ phi - u[3]

    u, _, ier, msg = fsolve(butler_res, np.concatenate((np.log(xs[1:]), [gam])), full_output=True, xtol=1e-14)
    if np.max(np.abs(butler_res(u))) > 1e-9:
        # fallback: longer, more damped fixed point, then Newton again
        for _ in range(600):
            gam = brentq(lambda g: xs_given(g, gs).sum() - 1.0, -50.0, 50.0 + abs(dGc), xtol=1e-15)
            xs = xs_given(gam, gs)
            gs = 0.8 * gs + 0.2 * BETA * partial_excess(xs / xs.sum(), T)
            pbar = 0.8 * pbar + 0.2 * float(xs @ phi) / xs.sum()
        u, _, ier, msg = fsolve(butler_res, np.concatenate((np.log(xs[1:]), [gam])), full_output=True, xtol=1e-14)
    if np.max(np.abs(butler_res(u))) > 1e-8:
        raise RuntimeError(f'Butler did not converge at T={T:.3f} K, dGc={dGc:.4g}: {msg}')
    xs = np.concatenate(([1 - np.exp(u[:3]).sum()], np.exp(u[:3])))
    return xs, u[3], xs


def load_sweep():
    rows = {}
    with open(SWEEP, newline='') as f:
        for r in csv.DictReader(f, delimiter=';'):
            if r['G_C'] in ('', 'nan'):
                continue
            rows.setdefault(r['surface'], []).append({k: float(r[k]) for k in ('grad', 'DT', 'gamma', 'G_C', 'r')})
    return rows


def main():
    xb = bulk_fractions()
    print('Bulk liquid X^B:', dict(zip(EL, np.round(xb, 6))))
    print('A_i [m2/mol]:', {e: f'{A[e]:.4e}' for e in EL})
    print('gamma_i (Turnbull) [J/m2]:', {e: round(GAMMA_I[e], 4) for e in EL})
    print('G_i^E,B at T_L [J/mol]:', dict(zip(EL, np.round(partial_excess(xb, T_L), 1))))

    xs0, gam0, u0 = butler(T_L, xb)
    print(f'\nStandard Butler at T_L: gamma = {gam0:.5f} J/m2, X^S =', dict(zip(EL, np.round(xs0, 6))))

    out = []
    for surf, pts in load_sweep().items():
        u_std, u_sh, u_rel = u0, u0, u0
        for p in pts:
            T = T_L - p['DT']
            xs_std, g_std, u_std = butler(T, xb, 0.0, u_std)
            xs_sh, g_sh, u_sh = butler(T, xb, p['G_C'], u_sh)
            # configurational term counted from equilibrium: (f'/f)(gamma - gamma_0), zero at the start of nucleation
            G_C_rel = p['G_C'] * (1 - GAMMA_0 / p['gamma'])
            xs_rel, g_rel, u_rel = butler(T, xb, G_C_rel, u_rel)
            out.append([surf, p['grad'], p['DT'], T, p['r'], p['gamma'], p['G_C'], gam0, g_std, g_sh,
                        *xs_std, *xs_sh, G_C_rel, g_rel, *xs_rel])
    head = ['surface', 'grad', 'DT', 'T', 'r', 'gamma_nuc', 'G_C', 'gamma_butler_TL', 'gamma_butler_std',
            'gamma_butler_shift'] + [f'XS_std_{e}' for e in EL] + [f'XS_shift_{e}' for e in EL] + \
           ['G_C_rel', 'gamma_butler_rel'] + [f'XS_rel_{e}' for e in EL]
    path = os.path.join(OUT, 'butler_AlSiMgFe.csv')
    with open(path, 'w', newline='') as f:
        w = csv.writer(f, delimiter=';')
        w.writerow(head)
        w.writerows([[r[0]] + [f'{v:.6e}' for v in r[1:]] for r in out])
    print('written', path)
    for surf in ('elastic', 'elasto-plastic'):
        sub = [r for r in out if r[0] == surf]
        for g in (sub[0][1], 7902.38, 1e5, sub[-1][1]):
            r = min(sub, key=lambda q: abs(log(q[1] / g)))
            print(f'{surf:15s} gradT={r[1]:10.4g} DT={r[2]:7.4f} gamma_nuc={r[5]:8.4f} '
                  f'Butler std={r[8]:.5f} shift={r[9]:8.4f}  XS_shift(Si,Mg,Fe)={r[15]:.4e},{r[16]:.4e},{r[17]:.4e}')


if __name__ == '__main__':
    main()
