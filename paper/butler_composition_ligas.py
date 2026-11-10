"""Surface composition of the FCC_A1 nucleus from Butler's formulation for the thermodynamic equilibrium shift
(Ferreira, Materials 2026, 19, 4676, Section 3.2, Eqs. 34a-c), with every energy taken from the nucleation solution
along the thermal gradient (option A), for the three alloys of the paper.

Butler is not used to compute the interfacial energy: the energies come from the nucleation solution with the elastic
and/or the elasto-plastic (stored_total) Gurtin-Murdoch surface,
    gamma = gamma_SL(grad T),   T = T_L - Delta T(grad T),   G_k^{E,C}/A_k = X_k^S Delta G_C(grad T),
and Butler's equations give the composition of the interface:
  - solutes SL_k (Eqs. 34b-c), unknowns X_k^S:
        gamma = gamma_k + RT/A_k ln(X_k^S/X_k^B) + (1/A_k)[G_k^{E,S}(T,X^S) - G_k^{E,B}(T,X^B)] - X_k^S Delta G_C
  - SV = Al, the solvent, fraction 1 - sum X_k (Eq. 34a), which closes the system with the effective interfacial
    energy of aluminium as the unknown; Al takes no configurational term (G^{E,C}_SV = 0, adopted):
        gamma_SV = gamma - RT/A_Al ln[(1 - sum X_k^S)/(1 - sum X_k^B)] - (1/A_Al)[G_Al^{E,S} - G_Al^{E,B}]
    (with G^{E,C}_SV = X_Al^S A_Al Delta G_C, test '(0)', gamma_SV is negative, since Delta G_C = -4 gamma_SL here).
  - sum X^S = 1 holds by construction: X_Al^S = 1 - sum X_k^S.
Assumption tests (python butler_composition_ligas.py --variant all): (0) the term on every component, (a) its opposite
sign (no solution above ~1e3 K/m), (c) Delta G_C from the chain-rule derivative of f(theta) (negligible term, the solutes
exhaust the interface above ~2.5e3 K/m); figures fig_assumption_<test>_butler_composition_<alloy>.png.

Liquid excess Gibbs energy: Redlich-Kister + Muggianu with the liquid parameters read from COST507-modified.tdb
(binaries and ternaries absent from the database are zero: Cu-Nb, Fe-Nb, Mg-Nb, Al-Cu-Fe, Al-Cu-Mg, Al-Cu-Nb, ...);
surface excess G^{E,S} = beta G^{E,B}(T, X^S), beta = 0.83 (Tanaka); pure-component gamma_k by Turnbull's rule and molar
interface areas A_k = 1.091 N_A^{1/3} V_k^{2/3}. The Turnbull value of gamma_Al is not used: it is the reference line.

Nucleation inputs:
  Al-0.8Si-0.6Mg-0.2Fe  outputs_AlSiMgFe/sweep_AlSiMgFe.csv (paper model; elastic and elasto-plastic)
  Al-3Cu-5Nb-0.1Fe      nucleacao/tablefull_paper_B_MRS_2026_GM_elastic.csv and ..._GM_EP_stored_total_refix_*.csv
  Al-3Cu-0.5Mg-0.2Fe    nucleacao/resultados/MR_2025_..._2026_paper/tablefull_..._stored_total_refix_*.csv
                        (elasto-plastic only: no elastic table of this alloy was kept)
Writes outputs_AlSiMgFe/butler_composition_<alloy>.csv and Fig_butler_composition_<alloy>.png.
Usage: python butler_composition_ligas.py [AlSiMgFe AlCuNbFe AlCuMgFe]
"""
import csv
import math
import os
import re
import sys
from itertools import combinations
from math import log

import numpy as np
from scipy.optimize import brentq, fsolve

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, 'outputs_AlSiMgFe')
TDB = os.path.join(HERE, 'COST507-modified.tdb')

R, N_A = 8.314462618, 6.02214076e23
F_AREA, BETA = 1.091, 0.83
#           M [g/mol]  V_solid [m^3/mol]  dH_m [J/mol]  alpha (Turnbull)
ELEM = {'AL': (26.9815, 10.00e-6, 10710.0, 0.45),
        'SI': (28.0855, 12.06e-6, 50210.0, 0.32),
        'MG': (24.305, 14.00e-6, 8480.0, 0.45),
        'FE': (55.845, 7.09e-6, 13810.0, 0.45),
        'CU': (63.546, 7.11e-6, 13260.0, 0.45),
        'NB': (92.906, 10.83e-6, 30000.0, 0.45)}

ALLOYS = {
    'AlSiMgFe': dict(title='Al-0.8Si-0.6Mg-0.2Fe', T_L=925.32, wt={'SI': 0.8, 'MG': 0.6, 'FE': 0.2},
                     sweeps={'paper': os.path.join(OUT, 'sweep_AlSiMgFe.csv')}),
    'AlCuNbFe': dict(title='Al-3Cu-5Nb-0.1Fe', T_L=936.58, wt={'CU': 3.0, 'NB': 5.0, 'FE': 0.1},
                     sweeps={'elastic': os.path.join(ROOT, 'nucleacao', 'tablefull_paper_B_MRS_2026_GM_elastic.csv'),
                             'elasto-plastic': os.path.join(ROOT, 'nucleacao',
                                 'tablefull_paper_B_MRS_2026_GM_EP_stored_total_refix_sy0.038572_H0.csv')}),
    'AlCuMgFe': dict(title='Al-3Cu-0.5Mg-0.2Fe', T_L=922.608, wt={'CU': 3.0, 'MG': 0.5, 'FE': 0.2},
                     sweeps={'elasto-plastic': os.path.join(ROOT, 'nucleacao', 'resultados',
                                 'MR_2025_dev_full_GLH_AlCuMgFe_signal_beta_2nd_Sim_2026_paper',
                                 'tablefull_paper_MR_2025_GM_EP_stored_total_refix_sy0.0422712_H0.csv')}),
}


# ------------------------------------------------------------------------------------- COST507 liquid
def read_liquid_parameters(path=TDB):
    """{('AL','CU'): {0: expr, 1: expr, ...}, ('AL','MG','SI'): {...}} from the first temperature interval."""
    txt = open(path, encoding='latin-1').read()
    pars = {}
    for m in re.finditer(r'PARAMETER\s+[GL]\(LIQUID,([A-Z,]+);(\d+)\)\s+[\d.]+\s+(.*?);', txt, re.S):
        els = tuple(m.group(1).split(','))
        if len(els) < 2:
            continue
        pars.setdefault(els, {})[int(m.group(2))] = ' '.join(m.group(3).split())
    return pars


_CODE = {}


def evaluate(expr, T):
    code = _CODE.get(expr)
    if code is None:
        code = _CODE[expr] = compile(expr.replace('LOG(', 'log(').replace('LN(', 'log('), '<tdb>', 'eval')
    return eval(code, {'__builtins__': {}}, {'T': T, 'log': math.log})


class Liquid:
    def __init__(self, els, pars):
        self.els = els
        self.bin, self.ter = [], []
        for c in combinations(range(len(els)), 2):
            key = tuple(sorted(els[i] for i in c))
            if key in pars:
                i, j = sorted(c, key=lambda q: els[q])
                self.bin.append((i, j, pars[key]))
        for c in combinations(range(len(els)), 3):
            key = tuple(sorted(els[i] for i in c))
            if key in pars:
                i, j, k = sorted(c, key=lambda q: els[q])
                self.ter.append((i, j, k, pars[key]))

    def g_excess(self, x, T):
        g = 0.0
        for i, j, Ls in self.bin:
            g += x[i] * x[j] * sum(evaluate(e, T) * (x[i] - x[j]) ** k for k, e in Ls.items())
        for i, j, k, Ls in self.ter:
            d = (1 - x[i] - x[j] - x[k]) / 3
            v = (x[i] + d, x[j] + d, x[k] + d)
            g += x[i] * x[j] * x[k] * sum(evaluate(Ls.get(n, '0'), T) * v[n] for n in range(3)) \
                if len(Ls) > 1 else x[i] * x[j] * x[k] * evaluate(Ls[0], T)
        return g

    def partial(self, x, T, h=1e-7):
        grad = np.empty(len(x))
        for k in range(len(x)):
            xp, xm = x.copy(), x.copy()
            xp[k] += h
            xm[k] -= h
            grad[k] = (self.g_excess(xp, T) - self.g_excess(xm, T)) / (2 * h)
        return self.g_excess(x, T) + grad - x @ grad


# --------------------------------------------------------------------------------------------- Butler
class NoSolution(RuntimeError):
    pass


class Alloy:
    def __init__(self, key, pars):
        cfg = ALLOYS[key]
        self.key, self.title, self.T_L = key, cfg['title'], cfg['T_L']
        self.els = ('AL',) + tuple(cfg['wt'])
        n = {e: cfg['wt'][e] / ELEM[e][0] for e in cfg['wt']}
        n['AL'] = (100 - sum(cfg['wt'].values())) / ELEM['AL'][0]
        tot = sum(n.values())
        self.xb = np.array([n[e] / tot for e in self.els])
        self.A = np.array([F_AREA * N_A ** (1 / 3) * ELEM[e][1] ** (2 / 3) for e in self.els])
        self.g0 = np.array([ELEM[e][3] * ELEM[e][2] / (N_A ** (1 / 3) * ELEM[e][1] ** (2 / 3)) for e in self.els])
        self.liq = Liquid(self.els, pars)
        self.sweeps = cfg['sweeps']

    def standard(self, T):
        """Standard Butler (no configurational term): gamma and X^S, all components, pure gamma_i of Turnbull."""
        xb, A, g0 = self.xb, self.A, self.g0
        gb = self.liq.partial(xb, T)
        gs = BETA * gb
        for _ in range(200):
            xk = lambda g: xb * np.exp(A * (g - g0 - (gs - gb) / A) / (R * T))
            gam = brentq(lambda g: xk(g).sum() - 1.0, -10.0, 10.0, xtol=1e-15)
            xs = xk(gam) / xk(gam).sum()
            gs = 0.5 * gs + 0.5 * BETA * self.liq.partial(xs, T)

        def res(u):
            xs = np.concatenate(([1 - np.exp(u[:-1]).sum()], np.exp(u[:-1])))
            if xs[0] <= 0:
                return np.full(len(u), 1e3)
            return g0 + R * T / A * np.log(xs / xb) + (BETA * self.liq.partial(xs, T) - gb) / A - u[-1]

        u = fsolve(res, np.concatenate((np.log(xs[1:]), [gam])), xtol=1e-14)
        assert np.max(np.abs(res(u))) < 1e-8, 'standard Butler did not converge'
        return u[-1], np.concatenate(([1 - np.exp(u[:-1]).sum()], np.exp(u[:-1])))

    def composition(self, T, gamma, dGc, guess=None, sgn=1.0, sv_term=True):
        """X^S of the solutes from Eqs. 34b-c with gamma and Delta G_C imposed; then gamma_SV of Al from Eq. 34a.
        sgn = +1: the configurational term enters as -X_k^S Delta G_C (Eq. 34 as written, -G^{E,C});
        sgn = -1: opposite sign (+G^{E,C}). sv_term = False: G^{E,C}_SV = 0 (Al takes no configurational term)."""
        xb, A, g0, liq = self.xb, self.A, self.g0, self.liq
        dGc = sgn * dGc
        sol = range(1, len(xb))
        gb = liq.partial(xb, T)

        def solute_x(k, gs_k):
            f = lambda lx: g0[k] + R * T / A[k] * (lx - log(xb[k])) + (gs_k - gb[k]) / A[k] - np.exp(lx) * dGc - gamma
            if dGc <= 0:                       # monotonic in X: a single root
                if f(0.0) <= 0:
                    return 1.0
                if f(-200.0) >= 0:
                    return 0.0
                return np.exp(brentq(f, -200.0, 0.0, xtol=1e-15))
            # dGc > 0 (variant a): f rises then falls; take the dilute root, none if the maximum is below zero
            grid = np.linspace(-200.0, 0.0, 8001)
            fv = f(grid)
            up = np.nonzero(fv >= 0)[0]
            if len(up) == 0:
                raise NoSolution(f'no root for {self.els[k]}: max of the equation {fv.max():.4f} J/m2 < 0')
            j = up[0]
            return 1.0 if j == 0 and fv[0] >= 0 else np.exp(brentq(f, grid[j - 1], grid[j], xtol=1e-15))

        x = xb.copy() if guess is None else guess.copy()
        for _ in range(200):
            gs = BETA * liq.partial(x, T)
            new = x.copy()
            for k in sol:
                new[k] = solute_x(k, gs[k])
            s = new[1:].sum()
            if s >= 1:
                new[1:] *= 0.999 / s
            new[0] = 1 - new[1:].sum()
            x = 0.5 * x + 0.5 * new

        def res(u):
            xs = np.concatenate(([1 - np.exp(u).sum()], np.exp(u)))
            if xs[0] <= 0:
                return np.full(len(u), 1e3)
            gs = BETA * liq.partial(xs, T)
            return np.array([g0[k] + R * T / A[k] * log(xs[k] / xb[k]) + (gs[k] - gb[k]) / A[k] - xs[k] * dGc - gamma
                             for k in sol])

        u, _, ier, msg = fsolve(res, np.log(x[1:]), full_output=True, xtol=1e-14)
        if np.max(np.abs(res(u))) > 1e-8:
            raise RuntimeError(f'{self.key}: composition did not converge at T={T:.3f} K, gamma={gamma:.4g}: {msg}')
        xs = np.concatenate(([1 - np.exp(u).sum()], np.exp(u)))
        gs = BETA * liq.partial(xs, T)
        g_sv = gamma - R * T / A[0] * log(xs[0] / xb[0]) - (gs[0] - gb[0]) / A[0] + (xs[0] * dGc if sv_term else 0.0)
        return xs, g_sv

    def load(self, chain_rule=False):
        """{surface: [(grad, DT, gamma, G_C), ...]} sorted by gradient.
        chain_rule: Delta G_C = (gamma/f) r (df/dtheta)(dtheta/dr), with df/dtheta = 3 sin^3(theta), instead of the
        model's dimensionless term gamma*dfthetadr/f (paper sweep only: the tables have no dtheta/dr)."""
        if chain_rule:
            data = {}
            with open(self.sweeps['paper'], newline='') as f:
                for r in csv.DictReader(f, delimiter=';'):
                    if r['G_C'] in ('', 'nan') or r['dtheta2_dr'] in ('', 'nan'):
                        continue
                    th, gam = float(r['theta']), float(r['gamma'])
                    fth = 2 - 3 * math.cos(th) + math.cos(th) ** 3
                    gc = gam / fth * float(r['r']) * 3 * math.sin(th) ** 3 * float(r['dtheta2_dr'])
                    data.setdefault(r['surface'], []).append((float(r['grad']), float(r['DT']), gam, gc))
            return {s: sorted(v) for s, v in data.items()}
        data = {}
        for surf, path in self.sweeps.items():
            with open(path, newline='') as f:
                rows = [{k.strip(): v for k, v in r.items()} for r in csv.DictReader(f, delimiter=';')]
            if surf == 'paper':
                for r in rows:
                    if r['G_C'] not in ('', 'nan'):
                        data.setdefault(r['surface'], []).append(
                            (float(r['grad']), float(r['DT']), float(r['gamma']), float(r['G_C'])))
            else:
                for r in rows:
                    if r['GC'] not in ('', 'nan') and float(r['gam_hom']) > 0:
                        data.setdefault(surf, []).append(
                            (float(r['dTdr']), float(r['DT']), float(r['gam_hom']), float(r['GC'])))
        return {s: sorted(v) for s, v in data.items()}


# --------------------------------------------------------------------------------------------- run
# Assumption tests on the configurational term (written as fig_assumption_<tag>_butler_composition_<alloy>.png):
VARIANTS = {
    # adopted (author's choice, 2026-11-09; formerly test 'sv_term_b'): the configurational term acts on the solutes,
    # Al (SV) takes none, so its effective interfacial energy stays positive and close to gamma_SL
    'base':        dict(sgn=1.0, sv_term=False, chain_rule=False,
                        text=r'Configurational term $-X_k^S\,\Delta G_C$ on the solutes, $G^{E,C}_{SV}=0$ for Al'),
    'all_terms_0': dict(sgn=1.0, sv_term=True, chain_rule=False,
                        text=r'(0) configurational term on every component, Al included: $\gamma_{SV}<0$'),
    'term_sign_a': dict(sgn=-1.0, sv_term=True, chain_rule=False,
                        text=r'(a) opposite sign of the configurational term: $+X_k^S\,\Delta G_C$'),
    'df_dr_c':     dict(sgn=1.0, sv_term=True, chain_rule=True,
                        text=r'(c) $\Delta G_C$ with the chain-rule derivative $r\,(\partial f/\partial\theta)(\partial\theta/\partial r)$'),
}


def run(key, pars, variant='base'):
    al = Alloy(key, pars)
    V = VARIANTS[variant]
    if V['chain_rule'] and 'paper' not in al.sweeps:
        print(f'\n== {al.title}: variant {variant} needs dtheta/dr, only in the paper sweep (Al-0.8Si-0.6Mg-0.2Fe) - skipped')
        return []
    g_std, xs_std = al.standard(al.T_L)
    print(f'\n== {al.title}  T_L = {al.T_L} K   variant {variant}: {V["text"]}')
    print('   X^B          :', dict(zip(al.els, np.round(al.xb, 6))))
    print('   gamma_i      :', dict(zip(al.els, np.round(al.g0, 4))))
    print('   liquid terms :', [tuple(al.els[i] for i in b[:2]) for b in al.liq.bin],
          [tuple(al.els[i] for i in t[:3]) for t in al.liq.ter])
    print(f'   standard Butler at T_L: gamma = {g_std:.5f} J/m2, X^S =', dict(zip(al.els, np.round(xs_std, 6))))
    out = []
    for surf, pts in al.load(V['chain_rule']).items():
        guess, first_fail = None, None
        for grad, DT, gam, gc in pts:
            T = al.T_L - DT
            try:
                xs, g_sv = al.composition(T, gam, gc, guess, V['sgn'], V['sv_term'])
                guess = xs
            except RuntimeError as e:          # NoSolution or no convergence: recorded as NaN
                xs, g_sv, guess = np.full(len(al.els), np.nan), np.nan, None
                if first_fail is None:
                    first_fail = (grad, gam, str(e).split(': ')[-1][:90])
            out.append([surf, grad, DT, T, gam, gc, g_sv, *xs])
        if first_fail:
            n_bad = sum(1 for r in out if r[0] == surf and np.isnan(r[6]))
            print(f'   {surf}: NO SOLUTION from gradT = {first_fail[0]:.4g} K/m (gamma_SL = {first_fail[1]:.4f} J/m2) on, '
                  f'{n_bad} of {len(pts)} points; {first_fail[2]}')
    tag = '' if variant == 'base' else f'assumption_{variant}_'
    path = os.path.join(OUT, f'butler_composition_{tag}{key}.csv')
    with open(path, 'w', newline='') as f:
        w = csv.writer(f, delimiter=';')
        w.writerow(['surface', 'grad', 'DT', 'T', 'gamma_nuc', 'G_C', 'gamma_SV'] + [f'XS_{e}' for e in al.els])
        w.writerows([[r[0]] + [f'{v:.6e}' for v in r[1:]] for r in out])
    for surf in dict.fromkeys(r[0] for r in out):
        sub = [r for r in out if r[0] == surf and not np.isnan(r[6])]
        if not sub:
            continue
        for g in (sub[0][1], 1e3, 1e4, 1e5, sub[-1][1]):
            r = min(sub, key=lambda q: abs(log(q[1] / g)))
            print(f'   {surf:15s} gradT={r[1]:10.4g} gamma={r[4]:8.4f} G_C={r[5]:10.3e} gamma_SV={r[6]:8.4f}  X^S ' +
                  ' '.join(f'{e}={v:.3e}' for e, v in zip(al.els, r[7:])))
    g_sv = [r[6] for r in out if not np.isnan(r[6])]
    if g_sv:
        print(f'   gamma_SV range (solved points): {min(g_sv):.4f} .. {max(g_sv):.4f} J/m2 '
              f'({"all positive" if min(g_sv) > 0 else "NEGATIVE values"})')
    figure(al, out, xs_std, variant)
    return out


def paper_figures(pars, keys, first=13):
    """Paper versions of the adopted figures (no suptitle), redrawn from the CSVs: Fig<n>_butler_composition_<alloy>.png."""
    for n, key in enumerate(keys, start=first):
        al = Alloy(key, pars)
        _, xs_std = al.standard(al.T_L)
        with open(os.path.join(OUT, f'butler_composition_{key}.csv'), newline='') as f:
            out = [[r['surface']] + [float(r[k]) for k in list(r)[1:]] for r in csv.DictReader(f, delimiter=';')]
        figure(al, out, xs_std, 'base', name=f'Fig{n}_butler_composition_{key}.png', title=False)


def figure(al, out, xs_std, variant='base', name=None, title=True):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'serif', 'font.size': 11, 'mathtext.fontset': 'dejavuserif'})
    col = {'AL': '#2a7ae2', 'SI': '#eb6834', 'MG': '#1a9e5c', 'FE': '#8a3ffc', 'CU': '#c0392b', 'NB': '#7f6000'}
    nm = {'AL': 'Al (SV)', 'SI': 'Si', 'MG': 'Mg', 'FE': 'Fe', 'CU': 'Cu', 'NB': 'Nb'}
    style = {'elastic': ('-', 2.0), 'elasto-plastic': ('--', 1.6)}
    surfs = list(dict.fromkeys(r[0] for r in out))
    fig, ax = plt.subplots(1, 2, figsize=(13.2, 5.0))
    for surf in surfs:
        ls, lw = style[surf]
        sub = [r for r in out if r[0] == surf]
        g = [r[1] for r in sub]
        for k, e in enumerate(al.els):
            ax[0].plot(g, [r[7 + k] for r in sub], ls, color=col[e], lw=lw,
                       label=nm[e] if surf == surfs[0] else None)
        ax[1].plot(g, [r[6] for r in sub], ls, color='#2a7ae2', lw=lw, label=f'$\\gamma_{{SV}}$, {surf}')
    for k, e in enumerate(al.els):
        ax[0].axhline(xs_std[k], color=col[e], ls=':', lw=1)
    ax[1].axhline(al.g0[0], color='0.3', ls=':', lw=1.2, label="Turnbull $\\gamma_{Al}$")
    for surf in surfs:
        ax[0].plot([], [], 'k' + style[surf][0], label=surf)
    ax[0].plot([], [], 'k:', label='standard Butler ($T_L$)')
    ax[0].set_yscale('log')
    ax[0].set_ylabel('Interface composition $X_i^S$ [-]')
    ax[1].set_ylabel(r'Effective interfacial energy of Al, $\gamma_{SV}$ [J m$^{-2}$]')
    for a, t in zip(ax, ('(a)', '(b)')):
        a.set_xscale('log')
        a.set_xlabel(r'Thermal gradient $\nabla T$ [K m$^{-1}$]')
        a.grid(True, which='both', ls='--', lw=0.4, alpha=0.6)
        a.text(0.03, 0.95, t, transform=a.transAxes, fontweight='bold', va='top')
        a.legend(fontsize=9, loc='best')
    if title:
        fig.suptitle(f'{al.title} (FCC_A1): surface composition from Butler with the equilibrium shift, '
                     'energies from the nucleation solution\n' + VARIANTS[variant]['text'], fontsize=11)
    else:
        ax[0].set_title(al.title, fontsize=11)
        ax[1].set_title(al.title, fontsize=11)
    fig.tight_layout()
    name = name or (f'Fig_butler_composition_{al.key}.png' if variant == 'base'
                    else f'fig_assumption_{variant}_butler_composition_{al.key}.png')
    path = os.path.join(HERE, name)
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print('   written', path)


if __name__ == '__main__':
    # usage: python butler_composition_ligas.py [--variant base|term_sign_a|sv_term_b|df_dr_c|all] [alloys...]
    args = sys.argv[1:]
    if args[:1] == ['--paper-figures']:      # Figures 13-15 of the paper, from the CSVs of the adopted solution
        paper_figures(read_liquid_parameters(), ['AlCuNbFe', 'AlSiMgFe', 'AlCuMgFe'])
        sys.exit()
    variants = ['base']
    if '--variant' in args:
        i = args.index('--variant')
        v = args[i + 1]
        variants = list(VARIANTS) if v == 'all' else [v]
        del args[i:i + 2]
    pars = read_liquid_parameters()
    for variant in variants:
        for key in (args or list(ALLOYS)):
            run(key, pars, variant)
