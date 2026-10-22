"""
Analytical x numerical verification of the elasto-plastic Gurtin-Murdoch (GM) surface stress tensor
implemented in Nucleation_model_..._ElastoPlastic_GM.py.

The functions under test are taken from the model file itself (ast extraction of the function
definitions only, so the script body of the model is not executed).

Part 1 - constitutive level, tangent plane, principal frame:
    excess stress  S~ = diag(m + q, m - q),  m = C_A (E_A - E_A^p),  q = C_D (e_d - e_d^p)
    E_A = e1 + e2 (areal strain), e_d = (e1 - e2)/2 (deviatoric strain)
    T1 equibiaxial (sphere), T2 pure shear, T3 uniaxial in-plane strain, T4 cyclic equibiaxial
    analytical (closed form, piecewise) x numerical return mapping (N steps) x RK45 of the rate form
    + consistency of the analytical dgamma/dr used by the 1st-order Gibbs-Thomson tensor
Part 2 - radial boundary-value problem: elastic core (lambda, mu) + elasto-plastic GM surface,
    generalized Young-Laplace  sigma_rr(R) = -2 tau / R  (liquid at zero pressure),
    surface areal strain  E_A = 2 ln(r/r_eq) + 2 u(R)/R,  tau = -sigma_0 + s
    analytical (u = A rho, scalar equation for A) x 1D finite volumes + Newton + return mapping
"""
import ast
import os
import numpy as np
from math import log, exp, cos, pi
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib.pyplot as plt

MODEL = 'Nucleation_model_FCC_A1_new_paper_B_MRS_Meeting_2026_dev_full_GLH_AlCuNbFe_signal_beta_2nd_Sim_2026_paper_ElastoPlastic_GM.py'

# ---------------------------------------------------------------- functions under test
src = open(MODEL, encoding='utf-8').read()
tree = ast.parse(src)
wanted = {'gm_surface_return_map', 's_closed_form_sphere', 'areal_strain_GM', 'dEA_dr_GM',
          'new_surface_state', 'gamma_func_r', 'gm_ep_surface'}
ns = {'log': log, 'exp': exp, 'cos': cos, 'pi': pi}
for node in tree.body:
    if isinstance(node, ast.FunctionDef) and node.name in wanted:
        exec(compile(ast.Module(body=[node], type_ignores=[]), MODEL, 'exec'), ns)
assert wanted <= set(ns), wanted - set(ns)

# model constants (same values as in the model file)
sigma_0, gamma_0 = 1.3761, 0.296
lambda0, mu0 = -5.0, 1.6
lambda_, mu = 5.034e10, 2.594e10
req = 6.3e-6
C_A, C_D = sigma_0, 2.0 * mu0
ns.update(sigma_0=sigma_0, gamma_0=gamma_0, lambda0=lambda0, mu0=mu0, lambda_=lambda_, mu=mu,
          C_A=C_A, C_D=C_D, SURFACE_MODEL='elastoplastic', DISSIPATION_MODE='stored',
          PLASTIC_REFERENCE_UPDATE=False, sigma_yA0=0.5, H_A=0.1 * sigma_0, sigma_yD0=0.5, H_D=0.0)
return_map = ns['gm_surface_return_map']
closed_sphere = ns['s_closed_form_sphere']


# ---------------------------------------------------------------- analytical 1D linear isotropic hardening
def analytic_1d(eps_path, C, sy0, H):
    """Exact piecewise solution of 1D linear-elastic / linear isotropic hardening plasticity along a
    piecewise-linear strain path (vertices eps_path). Returns a function eps(t) -> stress on a fine grid."""
    out_e, out_s = [eps_path[0]], [0.0]
    s, ep, al = 0.0, 0.0, 0.0
    for a, b in zip(eps_path[:-1], eps_path[1:]):
        d = np.sign(b - a)
        # elastic until |s| reaches the current yield stress in the loading direction
        sy = sy0 + H * al
        e_yield = ep + d * sy / C            # strain at which s = d * sy
        if d * (e_yield - a) > 0 and d * (b - e_yield) > 0:
            # elastic part a -> e_yield, then plastic part e_yield -> b
            for e in np.linspace(a, e_yield, 50)[1:]:
                out_e.append(e); out_s.append(C * (e - ep))
            for e in np.linspace(e_yield, b, 200)[1:]:
                dal = abs(e - e_yield) * C / (C + H)
                out_e.append(e); out_s.append(d * (sy + H * dal))
            al += abs(b - e_yield) * C / (C + H)
            s = d * (sy0 + H * al)
            ep = b - s / C
        else:
            for e in np.linspace(a, b, 100)[1:]:
                out_e.append(e); out_s.append(C * (e - ep))
            s = C * (b - ep)
    return np.array(out_e), np.array(out_s)


def numeric_return_map(e1, e2, sy_A0, H_A, sy_D0, H_D):
    """Return mapping of the model along the principal-strain path (e1[k], e2[k])."""
    st = ns['new_surface_state']()
    m, q = [], []
    for a, b in zip(e1, e2):
        s_, q_, *_ = return_map(a + b, 0.5 * (a - b), st, C_A, C_D, sy_A0, H_A, sy_D0, H_D, plastic=True)
        m.append(s_); q.append(q_)
    return np.array(m), np.array(q)


def rate_form_rk45(E_of_t, t_end, C, sy0, H):
    """Independent numerical integration of the rate form ds/dt = k_t(s, alpha, dE/dt) dE/dt (RK45)."""
    def rhs(t, y):
        s, al = y
        dE = (E_of_t(t + 1e-9) - E_of_t(t - 1e-9)) / 2e-9
        sy = sy0 + H * al
        loading = abs(s) >= sy * (1 - 1e-9) and np.sign(dE) == np.sign(s)
        if loading:
            kt = C * H / (C + H)
            return [kt * dE, abs(dE) * C / (C + H)]
        return [C * dE, 0.0]
    sol = solve_ivp(rhs, (0, t_end), [0.0, 0.0], method='RK45', max_step=t_end / 4000, rtol=1e-10, atol=1e-13,
                    dense_output=True)
    return sol


print('=' * 92)
print('PART 1 - tensor in the tangent plane: analytical x numerical')
print('=' * 92)
sy_A0, H_A, sy_D0, H_D = 0.5, 0.1 * sigma_0, 0.5, 0.05 * C_D
report = []


def check(name, comp, e_num, s_num, e_an, s_an):
    s_ref = np.interp(e_num, e_an, s_an) if np.all(np.diff(e_an) > 0) else None
    return s_ref


# T1 equibiaxial (sphere): E_A = 2 eps, e_d = 0
for N in (10, 100, 1000):
    eps = np.linspace(0.0, -0.6, N + 1)
    m, q = numeric_return_map(eps, eps, sy_A0, H_A, sy_D0, H_D)
    m_an = np.array([closed_sphere(2 * e, C_A, sy_A0, H_A) for e in eps])
    report.append(('T1 equibiaxial', 'm (areal)', N, np.max(np.abs(m - m_an)), np.max(np.abs(q))))
# T2 pure shear: E_A = 0, e_d = g
for N in (10, 100, 1000):
    g = np.linspace(0.0, 0.6, N + 1)
    m, q = numeric_return_map(g, -g, sy_A0, H_A, sy_D0, H_D)
    q_an = np.array([closed_sphere(e, C_D, sy_D0, H_D) for e in g])
    report.append(('T2 pure shear', 'q (J2)', N, np.max(np.abs(q - q_an)), np.max(np.abs(m))))
# T3 uniaxial in-plane strain: E_A = eps, e_d = eps/2 (both mechanisms)
for N in (10, 100, 1000):
    eps = np.linspace(0.0, -0.8, N + 1)
    m, q = numeric_return_map(eps, 0 * eps, sy_A0, H_A, sy_D0, H_D)
    m_an = np.array([closed_sphere(e, C_A, sy_A0, H_A) for e in eps])
    q_an = np.array([closed_sphere(e / 2, C_D, sy_D0, H_D) for e in eps])
    report.append(('T3 uniaxial', 'm and q', N, max(np.max(np.abs(m - m_an)), np.max(np.abs(q - q_an))), float('nan')))
# T4 cyclic equibiaxial: E_A 0 -> -1.2 -> +0.6 -> -1.6
verts = [0.0, -1.2, 0.6, -1.6]
e_an, s_an = analytic_1d(verts, C_A, sy_A0, H_A)
for N in (30, 300, 3000):
    EA = np.concatenate([np.linspace(a, b, N + 1)[:-1] for a, b in zip(verts[:-1], verts[1:])] + [[verts[-1]]])
    m, q = numeric_return_map(EA / 2, EA / 2, sy_A0, H_A, sy_D0, H_D)
    # the analytical solution is linear between its vertices (yield points included), so it is
    # evaluated exactly at the return-mapping points by interpolation along the path length
    seg_an = np.concatenate([[0], np.cumsum(np.abs(np.diff(e_an)))])
    seg_nu = np.concatenate([[0], np.cumsum(np.abs(np.diff(EA)))])
    s_an_on_nu = np.interp(seg_nu, seg_an, s_an)
    report.append(('T4 cyclic', 'm (areal)', N, np.max(np.abs(m - s_an_on_nu)), float('nan')))

print(f'{"test":<16}{"component":<12}{"steps":>7}{"max |num - analytical| [N/m]":>32}{"inactive comp.":>16}')
for name, comp, N, err, other in report:
    print(f'{name:<16}{comp:<12}{N:>7}{err:>32.3e}{other:>16.2e}')

# RK45 of the rate form for T4 (independent numerical method)
seg = np.abs(np.diff(verts)); t_k = np.concatenate([[0], np.cumsum(seg)])
E_of_t = lambda t: float(np.interp(t, t_k, verts))
sol = rate_form_rk45(E_of_t, t_k[-1], C_A, sy_A0, H_A)
seg_an = np.concatenate([[0], np.cumsum(np.abs(np.diff(e_an)))])
s_rk = sol.sol(seg_an)[0]
err_rk = np.max(np.abs(s_rk - s_an))
print(f'{"T4 cyclic":<16}{"m, RK45":<12}{"adapt.":>7}{err_rk:>32.3e}')

# consistency of the analytical dgamma/dr (input of Gamma 1st order) with gamma(r)
print('\ndgamma/dr used by Gamma^1st: analytical x central finite differences of gamma(r) and x integral')
for label, upd, sy, H, mode in [('elastic', None, None, None, 'stored'),
                                ('EP fixed ref, stored', False, 0.5, 0.1 * sigma_0, 'stored'),
                                ('EP fixed ref, thermo', False, 0.5, 0.1 * sigma_0, 'thermo'),
                                ('EP fixed ref, barrier', False, 0.5, 0.1 * sigma_0, 'barrier'),
                                ('EP updated ref, stored', True, 0.0695, 0.0, 'stored'),
                                ('EP updated ref, thermo', True, 0.05, 0.01 * sigma_0, 'thermo')]:
    ns['PLASTIC_REFERENCE_UPDATE'] = bool(upd)
    model = 'elastic' if upd is None else 'elastoplastic'
    pp = None if upd is None else (sy, H, sy_D0, 0.0)
    r = np.linspace(req, 0.35 * req, 6501)
    st = ns['new_surface_state']()
    g, dg = [], []
    for ri in r:
        _, _, gi, dgi, _, _ = ns['gm_ep_surface'](req, ri, st, model=model, mode=mode, plastic_params=pp)
        g.append(gi); dg.append(dgi)
    if upd is not None:
        bal = abs(st['wp'] - st['psi_p'] - st['q_diss']) / max(st['wp'], 1e-30)
        print(f'  {label:<22} energy balance |w_p - psi_p - q_diss| / w_p = {bal:.1e}  '
              f'(w_p = {st["wp"]:.4f}, psi_p = {st["psi_p"]:.4f}, q_diss = {st["q_diss"]:.4f} J/m^2)')
    g, dg = np.array(g), np.array(dg)
    fd = np.gradient(g, r)
    smooth = np.ones_like(r, bool)
    if upd is not None:   # exclude the two points around the yield onset, where dgamma/dr is discontinuous
        k = np.argmax(np.abs(np.diff(dg)) > 0.2 * np.abs(dg[1:]))
        smooth[max(k - 1, 0):k + 3] = False
    # scaled by max|dgamma/dr|: after yield with H = 0 and the updated reference dgamma/dr -> 0
    rel_fd = np.max(np.abs(fd[1:-1][smooth[1:-1]] - dg[1:-1][smooth[1:-1]])) / np.max(np.abs(dg))
    g_int = g[0] + np.concatenate([[0], np.cumsum(0.5 * (dg[1:] + dg[:-1]) * np.diff(r))])
    rel_int = np.max(np.abs(g_int - g) / g)
    print(f'  {label:<22} max |FD - analytical| / max|dgamma/dr| = {rel_fd:.2e}   '
          f'max rel |gamma - (gamma_eq + int dgamma)| = {rel_int:.2e}')
ns['PLASTIC_REFERENCE_UPDATE'] = False

# ---------------------------------------------------------------- PART 2
print('\n' + '=' * 92)
print('PART 2 - radial BVP: elastic core + elasto-plastic GM surface + Young-Laplace')
print('=' * 92)
K3 = 3 * lambda_ + 2 * mu
SY, HA = 0.5, 0.1 * sigma_0


def surface_closed(EA):
    return closed_sphere(EA, C_A, SY, HA)


def bvp_analytic(r):
    """u = A rho in the core; sigma_rr = (3 lambda + 2 mu) A = -2 tau / R; tau = -sigma_0 + s(E_A)."""
    Egeo = 2 * log(r / req)
    F = lambda A: K3 * A + 2 * (-sigma_0 + surface_closed(Egeo + 2 * A)) / r
    A = brentq(F, -1e-3, 1e-3, xtol=1e-22, rtol=1e-15)
    s = surface_closed(Egeo + 2 * A)
    return A, s, -sigma_0 + s


def bvp_fv(r_path, N):
    """1D finite volumes in rho on [0, R], Newton on (u_1..u_N), return mapping at the surface node.
    Surface plastic state is carried along the r-path (r decreasing); the core reference is the current R."""
    st = ns['new_surface_state']()
    res = []
    for R in r_path:
        rho = np.linspace(0, R, N + 1); h = R / N
        Egeo = 2 * log(R / req)
        u = np.zeros(N + 1)
        for it in range(30):
            trial = dict(st)
            EA = Egeo + 2 * u[-1] / R
            s, _, kA, *_ = return_map(EA, 0.0, trial, C_A, C_D, SY, HA, sy_D0, 0.0, plastic=True)
            tau = -sigma_0 + s
            # residual and Jacobian (linear elastic core; nonlinear only through tau at the surface)
            Rv = np.zeros(N); J = np.zeros((N, N))
            def srr_face(j):   # sigma_rr at face j+1/2 from u_j, u_{j+1}
                rf = 0.5 * (rho[j] + rho[j + 1])
                return rf, (lambda_ + 2 * mu) / h, lambda_ / rf
            for j in range(1, N + 1):
                # face j+1/2 (interior) or the surface traction (j = N)
                rl, a_l, b_l = srr_face(j - 1)
                sl = a_l * (u[j] - u[j - 1]) + b_l * (u[j] + u[j - 1])
                vol = h if j < N else h / 2
                rc = rho[j] if j < N else R - h / 4
                du_c = ((u[j + 1] - u[j - 1]) / (2 * h)) if j < N else (u[j] - u[j - 1]) / h
                stt = lambda_ * du_c + 2 * (lambda_ + mu) * (u[j] / rho[j])
                if j < N:
                    rr, a_r, b_r = srr_face(j)
                    sr = a_r * (u[j + 1] - u[j]) + b_r * (u[j + 1] + u[j])
                    Rv[j - 1] = rr ** 2 * sr - rl ** 2 * sl - 2 * rc * stt * vol
                else:
                    Rv[j - 1] = R ** 2 * (-2 * tau / R) - rl ** 2 * sl - 2 * rc * stt * vol
            # numerical Jacobian of the (small, banded) system, analytical surface tangent kA inside tau
            eps_u = 1e-15 + 1e-7 * max(1e-14, np.max(np.abs(u)))
            base = Rv.copy()
            for k in range(1, N + 1):
                u2 = u.copy(); u2[k] += eps_u
                trial2 = dict(st)
                s2, *_ = return_map(Egeo + 2 * u2[-1] / R, 0.0, trial2, C_A, C_D, SY, HA, sy_D0, 0.0, plastic=True)
                tau2 = -sigma_0 + s2
                R2 = np.zeros(N)
                for j in range(1, N + 1):
                    rl, a_l, b_l = srr_face(j - 1)
                    sl = a_l * (u2[j] - u2[j - 1]) + b_l * (u2[j] + u2[j - 1])
                    vol = h if j < N else h / 2
                    rc = rho[j] if j < N else R - h / 4
                    du_c = ((u2[j + 1] - u2[j - 1]) / (2 * h)) if j < N else (u2[j] - u2[j - 1]) / h
                    stt = lambda_ * du_c + 2 * (lambda_ + mu) * (u2[j] / rho[j])
                    if j < N:
                        rr, a_r, b_r = srr_face(j)
                        sr = a_r * (u2[j + 1] - u2[j]) + b_r * (u2[j + 1] + u2[j])
                        R2[j - 1] = rr ** 2 * sr - rl ** 2 * sl - 2 * rc * stt * vol
                    else:
                        R2[j - 1] = R ** 2 * (-2 * tau2 / R) - rl ** 2 * sl - 2 * rc * stt * vol
                J[:, k - 1] = (R2 - base) / eps_u
            du = np.linalg.solve(J, -base)
            u[1:] += du
            if np.max(np.abs(du)) < 1e-12 * max(1e-20, np.max(np.abs(u))) + 1e-30:
                break
        # commit the surface state at the converged displacement
        s, *_ = return_map(Egeo + 2 * u[-1] / R, 0.0, st, C_A, C_D, SY, HA, sy_D0, 0.0, plastic=True)
        res.append((R, u[-1] / R, s, -sigma_0 + s, (lambda_ + 2 * mu) * (u[1] - u[0]) / h + 2 * lambda_ * u[1] / rho[1],
                    np.max(np.abs(u[1:] / rho[1:] - u[-1] / R)) / abs(u[-1] / R), it + 1))
    return np.array(res)


r_path = req * np.linspace(1.0, 0.35, 66)[1:]
an = np.array([bvp_analytic(r) for r in r_path])
print(f'{"cells":>6}{"max rel err A=u(R)/R":>24}{"max |err s| [N/m]":>20}{"max |err tau| [N/m]":>21}{"max nonlin. of u":>18}{"Newton it.":>12}')
fv_runs = {}
for N in (8, 16, 32, 64):
    fv = bvp_fv(r_path, N)
    fv_runs[N] = fv
    errA = np.max(np.abs(fv[:, 1] - an[:, 0]) / np.abs(an[:, 0]))
    errs = np.max(np.abs(fv[:, 2] - an[:, 1]))
    errt = np.max(np.abs(fv[:, 3] - an[:, 2]))
    print(f'{N:>6}{errA:>24.3e}{errs:>20.3e}{errt:>21.3e}{np.max(fv[:, 5]):>18.2e}{int(np.max(fv[:, 6])):>12}')

# bulk partition: BVP surface strain x the model's E_A = 2 ln(y/y0)
EA_model = np.array([ns['areal_strain_GM'](req, r, lambda_, mu, lambda0, mu0) for r in r_path])
EA_bvp = 2 * np.log(r_path / req) + 2 * an[:, 0]
print(f'\nbulk partition of the areal strain: |E_A(BVP) - E_A(model, 2 ln y/y0)| <= {np.max(np.abs(EA_bvp - EA_model)):.2e} '
      f'(|E_A| up to {np.max(np.abs(EA_model)):.2f}); core strain A = u(R)/R up to {np.max(np.abs(an[:, 0])):.2e}')

# ---------------------------------------------------------------- figure
fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
eps = np.linspace(0.0, -0.8, 41)
m, q = numeric_return_map(eps, 0 * eps, sy_A0, H_A, sy_D0, H_D)
e_f = np.linspace(0, -0.8, 400)
ax = axs[0]
ax.plot(-e_f, [-closed_sphere(e, C_A, sy_A0, H_A) for e in e_f], color='#2a78d6', lw=2, label='analytical, mean $-m$')
ax.plot(-e_f, [-closed_sphere(e / 2, C_D, sy_D0, H_D) for e in e_f], color='#eb6834', lw=2, ls='--', label='analytical, deviator $-q$')
ax.plot(-eps, -m, 'o', color='#2a78d6', ms=7, mfc='none', label='return mapping (40 steps)')
ax.plot(-eps, -q, 's', color='#eb6834', ms=7, mfc='none', label='return mapping, deviator')
ax.set_xlabel(r'compressive in-plane strain $-\varepsilon_1$  ($\varepsilon_2 = 0$)', fontsize=11)
ax.set_ylabel(r'excess surface stress [N m$^{-1}$]', fontsize=11)
ax.set_title('(a) T3 uniaxial in-plane strain: both mechanisms', fontsize=11, loc='left')
ax.legend(fontsize=9, frameon=True, edgecolor='k')
ax = axs[1]
ax.plot(seg_an, s_an, color='#2a78d6', lw=2, label='analytical (piecewise)')
ax.plot(seg_an[::12], s_rk[::12], '^', color='#1baf7a', ms=6, mfc='none', label='RK45 of the rate form')
N = 30
EA = np.concatenate([np.linspace(a, b, N + 1)[:-1] for a, b in zip(verts[:-1], verts[1:])] + [[verts[-1]]])
mm, _ = numeric_return_map(EA / 2, EA / 2, sy_A0, H_A, sy_D0, H_D)
ax.plot(np.concatenate([[0], np.cumsum(np.abs(np.diff(EA)))]), mm, 'o', color='#eb6834', ms=5, mfc='none', label='return mapping (30 steps/branch)')
ax.set_xlabel(r'accumulated path length $\int |dE_A|$', fontsize=11)
ax.set_ylabel(r'mean excess stress $m$ [N m$^{-1}$]', fontsize=11)
ax.set_title(r'(b) T4 cyclic equibiaxial: $E_A$ 0 $\to$ -1.2 $\to$ 0.6 $\to$ -1.6', fontsize=11, loc='left')
ax.legend(fontsize=9, frameon=True, edgecolor='k')
ax = axs[2]
ax.plot(1 - r_path / req, an[:, 1], color='#2a78d6', lw=2, label=r'analytical, $s$')
ax.plot(1 - r_path / req, fv_runs[32][:, 2], 'o', color='#2a78d6', ms=5, mfc='none', label=r'finite volumes (32 cells), $s$')
ax.plot(1 - r_path / req, an[:, 2], color='#eb6834', lw=2, ls='--', label=r'analytical, $\tau = -\sigma_0 + s$')
ax.plot(1 - r_path / req, fv_runs[32][:, 3], 's', color='#eb6834', ms=5, mfc='none', label=r'finite volumes, $\tau$')
ax.set_xlabel(r'$1 - r/r_{eq}$', fontsize=11)
ax.set_ylabel(r'surface stress [N m$^{-1}$]', fontsize=11)
ax.set_title('(c) radial BVP: core + EP surface + Young-Laplace', fontsize=11, loc='left')
ax.legend(fontsize=9, frameon=True, edgecolor='k')
for a in axs:
    a.grid(True, linestyle='--', color='darkgray', linewidth=0.5)
fig.suptitle(r'Elasto-plastic Gurtin-Murdoch surface stress: analytical x numerical ($\sigma_y^s$ = 0.5 N/m, $H_A$ = 0.1$\sigma_0$, $\sigma_{yD}$ = 0.5 N/m, $H_D$ = 0.05$C_D$)', fontsize=12)
fig.tight_layout()
fig.savefig('Figure_GM_EP_tensor_analytical_vs_numerical.png', dpi=200, bbox_inches='tight')
print('\nfigure: Figure_GM_EP_tensor_analytical_vs_numerical.png')
