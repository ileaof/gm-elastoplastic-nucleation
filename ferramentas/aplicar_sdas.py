"""
Copies the SDAS (MRB) scripts into modelos_elasto-plastos/sdas and adds the elasto-plastic 2nd-order
tensor, keeping the elastic 1st- and 2nd-order tensors of the original. The original files are only read.

Changes in each copy:
  - Gibbs_Thomson_het_2nd_EP: 2nd-order tensor at the mean gradient from the elasto-plastic copy of the
    nucleation model of the same alloy (gamma_2nd_EP_ligas.txt, written by rodar_modelos_EP.py)
  - SDAS_EP_Eq = calc(Gibbs_Thomson_het_2nd_EP, ...), its power-law fit, a curve in Fig. 1 after the
    curves of the original and its label before the experimental label
  - the curves, labels and legend of the original are kept as they are
  - Al-3Cu-5Nb-0.1Fe (local G_L measured): a second elasto-plastic curve with the elastic restoration driven
    by the rate of density change, Gamma(P) from gamma_restaurado_AlCuNbFe.csv (ajustar_restauracao.py)
Usage:  python aplicar_sdas.py   (after rodar_modelos_EP.py)
"""
import os
import re
import shutil

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SRC = os.path.join(ROOT, 'Sim Luane Final')
DST = os.path.join(ROOT, 'modelos_elasto-plastos', 'sdas')
HERE = os.path.dirname(__file__)

# SDAS script -> (alloy, nucleation model whose elasto-plastic copy gives Gamma_2nd_EP)
SDAS = {
    'scriptsdasmodelAl3Cu5Nb01Fe_zoom_horizontal_B_MRS_Meeting_2026.py':
        ('AlCuNbFe', 'Nucleation_model_FCC_A1_new_paper_B_MRS_Meeting_2026_dev_full_GLH_AlCuNbFe_signal_beta_2nd_Sim_2026_paper.py'),
    'scriptsdasmodelAl08Si06Mg02Fe_zoom_horizontal_new_script.py':
        ('AlSiMgFe', 'Nucleation_model_FCC_A1_new_paper_2024_dev_full_GLH_AlSiMgFe_signal_beta_2nd_Sim_Luane.py'),
    'scriptsdasmodelAl08Si06Mg02Fe_zoom_horizontal_new_script_Ferreira_2024.py':
        ('AlSiMgFe', 'Nucleation_model_FCC_A1_new_paper_2024_dev_full_GLH_AlSiMgFe_signal_beta_2nd_Sim_Luane.py'),
    'scriptsdasmodelAl3Cu05Mg02Fe_zoom_horizontal_new_script_Ferreira_2025.py':
        ('AlCuMgFe', 'Nucleation_model_FCC_A1_new_paper_MR_2025_dev_full_GLH_AlCuMgFe_signal_beta_2nd_Sim_2026_paper.py'),
    'scriptsdasmodelAl3Cu05Mg02Fe_zoom_horizontal_new_script_Ferreira_2025_dev_gam.py':
        ('AlCuMgFe', 'Nucleation_model_FCC_A1_new_paper_MR_2025_dev_full_GLH_AlCuMgFe_signal_beta_2nd_Sim_dev_gam.py'),
}
HELPERS = ['midpoint.py']


def short_name(nucleation_file):
    tag = nucleation_file[:-3].replace(' ', '_')
    return tag.replace('Nucleation_model_FCC_A1_new_paper_', '').replace('Nucleation_model_', '')[:60]


def read_table(name):
    out = {}
    for line in open(os.path.join(HERE, name)):
        if line.strip() and not line.startswith('#'):
            p = line.split()
            out[p[0]] = p[1:]
    return out


def insert_after(pattern, code, text, what):
    m = re.search(pattern, text, flags=re.M)
    if not m:
        raise RuntimeError(f'pattern not found ({what}): {pattern}')
    return text[:m.end()] + '\n' + code + text[m.end():]


def sub_once(pattern, repl, text, what):
    new, k = re.subn(pattern, repl, text, count=1, flags=re.M)
    if k != 1:
        raise RuntimeError(f'pattern not found ({what}): {pattern}')
    return new


RESTORATION_CSV = 'gamma_restaurado_AlCuNbFe.csv'
# only the elastic part |s| <= sigma_y is restored; the permanent deformation stays stored (author, 2026-10-22)
RESTORATION_IN_SDAS = True


def read_restoration():
    out = {}
    for line in open(os.path.join(HERE, 'restauracao_AlCuNbFe.txt')):
        if line.strip() and not line.startswith('#'):
            k, v = line.split()[:2]
            out[k] = float(v)
    return out


def add_restoration(t):
    """Elasto-plastic 2nd-order tensor with the elastic restoration driven by the rate of density change."""
    rs = read_restoration()
    code = '\n'.join([
        '# --- elastic restoration driven by the rate of density change (ferramentas/ajustar_restauracao.py) ---',
        '# eps_v = beta_s T_dot/DT_SL, f = eps_v tau_r/(1 + eps_v tau_r) (f = 1 without tau_r), Gamma = Gamma_EP - f |s|/|DS_V|',
        f"# beta_s = {rs['beta_s']:.4g}, DT_SL = {rs['DT_SL']:.4g} K (mean integral 5-90 mm), tau_r = {rs['tau_r']:.4g} s",
        '# constant Gamma at the mean experimental gradient, as the other MRB curves, minus the restorable part |s|/|DS_V|',
        f"Gibbs_Thomson_het_2nd_EP_rest = Gibbs_Thomson_het_2nd_EP - {rs['dGamma_el_mean']:.6e}   # [m.K]",
        'SDAS_EP_rest_Eq, beta_EPr_1, beta_EPr_2, beta_EPr_3 = calc(Gibbs_Thomson_het_2nd_EP_rest, delta1m_Eq, delta2m_Eq, '
        'delta3m_Eq, kef1_Eq, kef2_Eq, kef3_Eq, w1_Eq, w2_Eq, w3_Eq)',
        'z3n_EP_rest = np.poly1d(np.polyfit(np.log(tSL), np.log(SDAS_EP_rest_Eq*1e06), 1))'])
    t = insert_after(r'^z3n_EP = np\.poly1d\(z3_EP\).*$', code, t, 'z3n_EP')
    t = insert_after(r"^ax1\.plot\(tSL,\(np\.exp\(z3n_EP\[0\]\)\)\*tSL\*\*z3n_EP\[1\],':k'.*$",
                     "ax1.plot(tSL,(np.exp(z3n_EP_rest[0]))*tSL**z3n_EP_rest[1],'-',color='0.45',linewidth = 2.5)   # elasto-plastic, 2nd order, elastic restoration",
                     t, 'elasto-plastic curve of Fig. 1')
    label = ("r'MRB Model with $\\Gamma ^{2nd-Order} _{Het}$ elasto-plastic, elastic part restored ($f$ = 1), "
             "$\\mathit{SDAS_{MRB}\\,=\\,%.2f\\,t_{SL}^{%.3f} \\,\\, \\mathrm{[\\mu m]}}$' % "
             "(np.exp(z3n_EP_rest[0]), z3n_EP_rest[1]), ")
    m = re.search(r"^plt\.legend\(\(r'RB SDAS.*$", t, flags=re.M)
    leg = m.group(0)
    k = leg.find("r'Experimental")
    return t[:m.start()] + leg[:k] + label + leg[k:] + t[m.end():]


def convert_safe(text, gt2, sy, copy_name, alloy=None):
    """Adds the elasto-plastic curve after the last MRB curve of Fig. 1; the original curves are untouched."""
    t = convert_header(text, gt2, sy, copy_name)
    t = sub_once(r"^(ax1\.plot\(tSL,\(np\.exp\(z3n_1st_Eq\[0\]\)\)\*tSL\*\*z3n_1st_Eq\[1\],'--k',linewidth = 1\.5\).*)$",
                 r'\1' "\nax1.plot(tSL,(np.exp(z3n_EP[0]))*tSL**z3n_EP[1],':k',linewidth = 2.5)   # elasto-plastic, 2nd order",
                 t, 'MRB curve of Fig. 1')
    label = ("r'MRB Model with $\\Gamma ^{2nd-Order} _{Het}$ elasto-plastic ($\\sigma _y ^s$ = %g N/m), "
             "$\\mathit{SDAS_{MRB}\\,=\\,%.2f\\,t_{SL}^{%.3f} \\,\\, \\mathrm{[\\mu m]}}$' % "
             f"({sy:g}, np.exp(z3n_EP[0]), z3n_EP[1]), ")
    m = re.search(r"^plt\.legend\(\(r'RB SDAS.*$", t, flags=re.M)
    if not m:
        raise RuntimeError('legend of Fig. 1 not found')
    leg = m.group(0)
    k = leg.find("r'Experimental")
    if k < 0:
        raise RuntimeError('experimental label not found in the legend')
    t = t[:m.start()] + leg[:k] + label + leg[k:] + t[m.end():]
    restored = RESTORATION_IN_SDAS and alloy == 'AlCuNbFe' and os.path.exists(os.path.join(HERE, RESTORATION_CSV))
    if restored:
        t = add_restoration(t)
    return legend_handles(t, restored)


def legend_handles(t, restored):
    """Fig. 1: each label of the legend next to the line of its own curve (the labels are matched to the lines in
    plotting order otherwise, and the 2nd-order curve is plotted before the 1st-order one)."""
    for name, pat in (('_h_RB', r"^ax1\.plot\(tSL,\(np\.exp\(z3n_RB\[0\]\)\)"),
                      ('_h_2nd', r"^ax1\.plot\(tSL,\(np\.exp\(z3n_Eq\[0\]\)\)"),
                      ('_h_1st', r"^ax1\.plot\(tSL,\(np\.exp\(z3n_1st_Eq\[0\]\)\)"),
                      ('_h_EP', r"^ax1\.plot\(tSL,\(np\.exp\(z3n_EP\[0\]\)\)"),
                      ('_h_exp', r"^ax1\.plot\(xsc, ysc, 'sk'")):
        t = sub_once('(' + pat + ')', name + r', = \1', t, name)
    hs = '[_h_RB, _h_1st, _h_2nd, _h_EP' + (', _h_R' if restored else '') + ', _h_exp]'
    if restored:
        t = sub_once(r"^(ax1\.plot\(tSL,\(np\.exp\(z3n_EP_rest\[0\]\)\))", r'_h_R, = \1', t, '_h_R')
    return sub_once(r"^plt\.legend\(\(r'RB SDAS", 'plt.legend(' + hs + r", (r'RB SDAS", t, 'legend handles')


def convert_header(text, gt2, sy, copy_name):
    t = sub_once(r'^(SDAS_het_1st_Eq, .*)$',
                 r'\1' '\n# --- elasto-plastic Gurtin-Murdoch surface (modelos_elasto-plastos) ---\n'
                 f'# 2nd-order tensor at the mean gradient from nucleacao/{copy_name}_EP.py\n'
                 f'# (sigma_y = {sy:g} N/m, plastic area updating the reference, H = 0, loaded surface)\n'
                 f'Gibbs_Thomson_het_2nd_EP = {gt2:.6e}   # [m.K]\n'
                 'SDAS_EP_Eq, beta_EP_1, beta_EP_2, beta_EP_3 = calc(Gibbs_Thomson_het_2nd_EP, delta1m_Eq, delta2m_Eq, '
                 'delta3m_Eq, kef1_Eq, kef2_Eq, kef3_Eq, w1_Eq, w2_Eq, w3_Eq)\n'
                 "print(f'SDAS_EP_Eq (Gamma 2nd order, elasto-plastic = {Gibbs_Thomson_het_2nd_EP:.4e} m.K)')",
                 text, 'SDAS_het_1st_Eq')
    t = sub_once(r'^(z3n_1st_Eq\s*=\s*np\.poly1d\(z3_1st_Eq\).*)$',
                 r'\1' '\nlx_EP = np.log(tSL)\nly_EP = np.log(SDAS_EP_Eq*1e06)\nz3_EP = np.polyfit(lx_EP, ly_EP, 1)\nz3n_EP = np.poly1d(z3_EP)',
                 t, 'z3n_1st_Eq')
    return t


def main():
    os.makedirs(DST, exist_ok=True)
    for h in HELPERS:
        shutil.copy2(os.path.join(SRC, h), DST)
    gam = read_table('gamma_2nd_EP_ligas.txt')
    sig = read_table('sigma_y_ligas.txt')
    for sdas, (alloy, nucl) in SDAS.items():
        short = short_name(nucl)
        gt2 = float(gam[short][2])
        sy = float(sig[alloy][0])
        text = open(os.path.join(SRC, sdas), encoding='utf-8').read()
        new = convert_safe(text, gt2, sy, nucl[:-3], alloy)
        header = ('"""\nELASTO-PLASTIC COPY (modelos_elasto-plastos). Original: Sim Luane Final/%s\n'
                  'Adds the 2nd-order tensor of the elasto-plastic Gurtin-Murdoch surface (%s, sigma_y = %g N/m)\n'
                  'to the elastic 1st- and 2nd-order tensors of the original. See ../README.md.\n"""\n' % (sdas, alloy, sy))
        name = sdas[:-3] + '_EP.py'
        open(os.path.join(DST, name), 'w', encoding='utf-8').write(header + new)
        print(f'{sdas}  ->  sdas/{name}  ({alloy}: Gamma_2nd_EP = {gt2:.4e} m K, sigma_y = {sy:g} N/m)')


if __name__ == '__main__':
    main()
