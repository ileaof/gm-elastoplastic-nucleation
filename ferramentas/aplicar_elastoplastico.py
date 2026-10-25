"""
Copies the elastic nucleation models into modelos_elasto-plastos/nucleacao and applies the
elasto-plastic Gurtin-Murdoch surface (gm_elastoplastic.py). The original files are only read.

Changes in each copy (every pattern must be found, otherwise the copy is not written):
  - import of gm_elastoplastic
  - arrays enlarged (n = 400) and the r-sweep extended until Delta T > GM_DT_MAX (default 5 K)
  - sigma_func / gamma_func_r / dgammadr_func_r of the main loop -> gm_surf.step(r)
  - Delta T in closed form (the root that brentq(f_DT, 1e-4, 1) returns in the original scripts),
    gradient dTdr = Delta T / (8 pi r)
  - 1st-order theta: NaN when it has no root; 2nd-order theta: the sweep stops when it has no root
  - index j of the experimental gradient found from the sweep (was a fixed number)
  - output CSV names tagged with the surface model
Usage:  python aplicar_elastoplastico.py
"""
import os
import re
import shutil

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SRC = os.path.join(ROOT, 'Sim Luane Final')
DST = os.path.join(ROOT, 'modelos_elasto-plastos', 'nucleacao')

# original file (relative to SRC) -> (copy name, sigma_y of the alloy [N/m])
MODELS = {
    'Nucleation_model_FCC_A1_new_paper_2024_dev_full_GLH_AlSiMgFe_signal_beta_2nd_Sim_Luane.py': 'AlSiMgFe',
    'Nucleation_model_FCC_A1_new_paper_B_MRS_Meeting_2026_dev_full_GLH_AlCuNbFe_signal_beta_2nd_Sim_2026_paper.py': 'AlCuNbFe',
    'Nucleation_model_FCC_A1_new_paper_MR_2025_dev_full_GLH_AlCuMgFe_signal_beta_2nd_Sim.py': 'AlCuMgFe',
    'Nucleation_model_FCC_A1_new_paper_MR_2025_dev_full_GLH_AlCuMgFe_signal_beta_2nd_Sim_2026_paper.py': 'AlCuMgFe',
    'Nucleation_model_FCC_A1_new_paper_MR_2025_dev_full_GLH_AlCuMgFe_signal_beta_2nd_Sim_dev.py': 'AlCuMgFe',
    'Nucleation_model_FCC_A1_new_paper_MR_2025_dev_full_GLH_AlCuMgFe_signal_beta_2nd_Sim_dev_gam.py': 'AlCuMgFe',
    'Nucleation_model_FCC_A1_new_paper_MR_2025_dev_full_GLH_AlCuMgFe_signal_beta_2nd_Sim_dev_old.py': 'AlCuMgFe',
    'Nucleation_model_FCC_A1_new_paper_MR_2025_dev_full_GLH_AlCuMgFe_signal_beta_2nd_Sim_old.py': 'AlCuMgFe',
    'Nucleation_model_H2O_vapour_solid_Sim_2026_paper.py': 'H2O',
    os.path.join('métrica', 'Nucleation_model_FCC_A1_new_paper_MR_2025_dev_full_GLH_AlCuMgFe_signal_beta_2nd_Sim - Copia.py'): 'AlCuMgFe',
    os.path.join('scientific notation', 'Nucleation_model_FCC_A1_new_paper_MR_2025_dev_full_GLH_AlCuMgFe_signal_beta_2nd_Sim_2026_paper scient notation.py'): 'AlCuMgFe',
}
NOT_CONVERTED = ['Nucleation_model_FCC_A1_new_paper_proposal_luana_GL_data_indev.py']
HELPERS = ['derivative_shim.py', 'midpoint.py', 'mpl_textguard.py']
# energy of the surface in the copies (author, 2026-10-22): elastic + permanent deformation, stored for nucleation
ENERGY_MODE = 'stored_total'   # permanent deformation stored, fixed reference, minus dissipated heat (author, 2026-10-22)

# areal yield stress of the interface for each alloy [N/m], fitted to its SDAS (see ../README.md)
SIGMA_Y = {'AlCuNbFe': 0.0695, 'AlSiMgFe': 0.0695, 'AlCuMgFe': 0.0695, 'H2O': 0.0695}
if os.path.exists(os.path.join(os.path.dirname(__file__), 'sigma_y_ligas.txt')):
    for line in open(os.path.join(os.path.dirname(__file__), 'sigma_y_ligas.txt')):
        if line.strip() and not line.startswith('#'):
            k, v = line.split()[:2]
            SIGMA_Y[k] = float(v)


# part of the alloy name that identifies the alloy of each copy in the texts of the figures
ALLOY_TEXT = {'AlCuNbFe': '5wt%Nb', 'AlSiMgFe': '0.6wt%Mg', 'AlCuMgFe': '0.5wt%Mg', 'H2O': None}

RATE_VS_GRADIENT = r'''
# --- nucleation rate against the thermal gradient (elasto-plastic copy; the rate against the undercooling is kept) ---
_np = __import__('numpy')
_G = dTdr_v[1:npointsm1]
_EP = gm_surf.model == 'elastoplastic'
_sfx = ', elasto-plastic' if _EP else ', elastic'
_ovG = globals().get('_ov')
if _ovG is None and _EP and not __import__('os').environ.get('GM_DUMP_RATE'):
    _ovG = elastic_rate(__file__)
figG, axG = plt.subplots(1)
axG.set_title(r'Nucleation Statistics against the Thermal Gradient, phase FCC_A1', {'color': 'black', 'fontsize': 13})
axG.plot(_G, stat_hom_2nd_v[1:npointsm1], 'k-', linewidth=1.5, label='Homogeneous' + _sfx)
axG.plot(_G, stat_het_2nd_v[1:npointsm1], 'k-.', linewidth=2.5, label='Heterogeneous' + _sfx)
if _EP and _ovG is not None:
    axG.plot(_ovG['G'], _ovG['Ihom'], color='#2a78d6', linewidth=1.5, label='Homogeneous, elastic')
    axG.plot(_ovG['G'], _ovG['Ihet'], color='#2a78d6', linewidth=2.0, linestyle='-.', label='Heterogeneous, elastic')
_Gy = None
if _EP and gm_surf.yielded_at:
    _Gy = _G[int(_np.argmin(_np.abs(rc_v[1:npointsm1] - gm_surf.yielded_at)))]
    axG.axvspan(_G[0], _Gy, color='#2a78d6', alpha=0.08, lw=0)
    axG.axvline(_Gy, color='k', linewidth=1.0, linestyle='-.', label=r'Yield (elastic range shaded), $\nabla T$ = %.0f K/m' % _Gy)
axG.axvline(dTdrj, color='0.4', linewidth=1.0, linestyle=':', label=r'Mean experimental $\overline{\nabla T}$ = %g K/m' % dTdrj)
axG.set_xscale('log')
axG.set_xlabel(r'Thermal Gradient, $\mathrm{\nabla T \,\, [K.m^{-1}]}$', {'color': 'black', 'fontsize': 16})
axG.set_ylabel(r'Nucleation Sites, $I \,\, \mathrm{[m^{-3}s^{-1}]}$', {'color': 'black', 'fontsize': 16})
plt.setp(axG.get_xticklabels(), fontsize=14)
plt.setp(axG.get_yticklabels(), fontsize=14)
axG.grid(True, which='both', linestyle='--', color='0.85', linewidth=0.5)
axG.legend(title='2nd-order statistics', fontsize=8.5, frameon=True, edgecolor='k', shadow=True, loc='upper right')
# --- Gibbs-Thomson, surface energy and surface tension against the thermal gradient, whole range of the sweep ---
_ovOK = _EP and _ovG is not None
for _var, _key, _yl, _tt, _lg, _u in (
        ('GT_het_2nd_v', 'GT', r'Gibbs-Thomson 2nd order, $\Gamma^{2nd}$ [m K]', 'Second-Order Gibbs-Thomson', True, 'm K'),
        ('gam_hom_v', 'gam', r'Surface energy, $\gamma_{SL}$ [J m$^{-2}$]', 'Surface Energy', True, 'J/m$^2$')):
    vs_gradient_figure(_G, _np.asarray(globals()[_var][1:npointsm1], float), _yl,
                       _tt + ' against the Thermal Gradient, phase FCC_A1', _EP,
                       y_el=_ovG[_key] if (_ovOK and _key in _ovG) else None, G_el=_ovG['G'] if _ovOK else None,
                       G_yield=_Gy, G_mean=dTdrj, log_y=_lg, unit=_u)
# surface tension composed of the elastic and plastic components:
#   tau = (sigma_0 + s^e) + |sigma_0| E_A^p + sigma_y   (elastic component capped at sigma_0 - sigma_y after yielding)
_h = _np.asarray(gm_surf.hist[1:npointsm1], float)
_ext = ((_G, _h[:, 0], r'Elastic component, $\sigma_0 + s^e$', dict(color='k', lw=1.2, ls=':')),
        (_G, _h[:, 1], r'Plastic component, $|\sigma_0|\,E_A^p$', dict(color='#8a3ffc', lw=1.2, ls='-.'))) if _EP else ()
vs_gradient_figure(_G, _h[:, 3], r'Surface tension, $\tau_{SL}$ [N m$^{-1}$]',
                   'Surface Tension against the Thermal Gradient, phase FCC_A1', _EP,
                   y_el=_ovG['tau'] if (_ovOK and 'tau' in _ovG) else None, G_el=_ovG['G'] if _ovOK else None,
                   G_yield=_Gy, G_mean=dTdrj, log_y=False, unit='N/m', extra=_ext,
                   main_label=(r'Total, elastic + plastic: $\sigma_0 + s^e + |\sigma_0|E_A^p + \sigma_y^s$' if _EP else None))
'''


def sub_once(pattern, repl, text, name, flags=re.M):
    new, k = re.subn(pattern, repl, text, count=1, flags=flags)
    if k != 1:
        raise RuntimeError(f'{name}: pattern not found: {pattern}')
    return new


def convert(text, alloy):
    sy = SIGMA_Y[alloy]
    t = text
    # surface energy form of the original gamma_func_r (the *_dev_gam script uses the metric form)
    gfun = re.search(r'^def gamma_func_r\(.*?(?=^def )', t, flags=re.M | re.S).group(0)
    form = 'metric_2sphere' if 'mconst_sph' in gfun else 'legacy'
    t = sub_once(r'^(from derivative_shim import derivative.*)$',
                 r'\1\nfrom gm_elastoplastic import GMElastoPlasticSurface, DT_capillary   # elasto-plastic GM surface', t, 'import')
    t = sub_once(r'^n = 100\b.*$', 'n = 400  # enlarged for the elasto-plastic sweep (was 100)', t, 'n')
    t = sub_once(r"^(file = open\()'([^']+)\.csv','w'\)",
                 '# elasto-plastic Gurtin-Murdoch surface (gm_elastoplastic.py), sigma_y of this alloy: %g N/m\n'
                 '# surface energy (%s): elastic + permanent surface deformation, energy stored for nucleation\n'
                 'gm_surf = GMElastoPlasticSurface(req, gamma_0, sigma_0, lambda_, mu, lambda0, mu0, sigma_y_default=%r, surface_form=%r,\n'
                 '                                 mode_default=%r)\n'
                 "print(f'Surface: {gm_surf.describe()}')\n"
                 r"\1f'\2_GM_{gm_surf.tag}.csv','w')" % (sy, ENERGY_MODE, sy, form, ENERGY_MODE), t, 'file')
    t = re.sub(r"^(filefull = open\()'([^']+)\.csv','w'\)", r"\1f'\2_GM_{gm_surf.tag}.csv','w')", t, count=1, flags=re.M)
    t = sub_once(r'^npoints = \d+.*$', 'npoints = n  # the sweep stops at Delta T > GM_DT_MAX (was a fixed number)', t, 'npoints')
    t = sub_once(r'^npointsm1 = npoints - 1.*\n', '', t, 'npointsm1')
    t = sub_once(r'^(\s+)r = req - i \* dreq.*$',
                 r'\g<0>\n\1if r <= 0.0:\n\1    npoints = i\n\1    break', t, 'r')
    t = sub_once(r'^(\s+)sigm, surfstress = sigma_func\(req, r, .*$',
                 r'\1sigm, surfstress, gam, dfgamdr_ana = gm_surf.step(r)   # elasto-plastic GM surface', t, 'sigma_func')
    t = sub_once(r'^(\s+)gam = gamma_func_r\(gamma_0, r, req, surfstress\).*$', r'\1# gam: from gm_surf.step(r)', t, 'gamma_func_r')
    t = sub_once(r'^(\s+)dfgamdr_ana = dgammadr_func_r\(req, r, .*$', r'\1# dfgamdr_ana: from gm_surf.step(r)', t, 'dgammadr_func_r')
    t = sub_once(r'^(\s+)dTdr = brentq\(f,1,\d+, .*$',
                 r'\1dTdr = DT_capillary(gam, Delta_Sv_hom, r) / (8*pi*r)   # root of f in the original script', t, 'dTdr')
    t = sub_once(r'^(\s+)DT = brentq\(f_DT,0\.0001,1, .*$',
                 r'\1DT = DT_capillary(gam, Delta_Sv_hom, r)   # root of f_DT in the original script\n'
                 r'\1if DT > gm_surf.DT_max:\n\1    npoints = i\n\1    break', t, 'DT')
    t = sub_once(r'^(\s+)(thet_2nd = brentq\(f_het_2nd,.*)$',
                 r'\1try:\n\1    \2\n\1except (ValueError, RuntimeError):\n\1    npoints = i   # no 2nd-order angle: end of the sweep\n\1    break',
                 t, 'thet_2nd')
    t = sub_once(r'^(\s+)(thet = brentq\(f_het,.*)$',
                 r"\1try:\n\1    \2\n\1except (ValueError, RuntimeError):\n\1    thet = float('nan')   # no 1st-order angle (dgamma/dr -> 0 after yield)",
                 t, 'thet')
    t = sub_once(r'^file\.close\(\)\s*$',
                 'file.close()\nnpointsm1 = npoints - 1\n'
                 "print(f'Sweep: {npoints} points, grad T up to {dTdr_v[npointsm1]:.1f} K/m, DT up to {DT_v[npointsm1]:.4f} K'"
                 " + (f', yield at r/req = {gm_surf.yielded_at/req:.4f}' if gm_surf.yielded_at else ''))",
                 t, 'file.close')
    t = sub_once(r'^j = \d+.*$', '# j: found below from the gradient sweep (was a fixed number)', t, 'j')
    t = sub_once(r'^jp1 = j \+ 1.*$', '', t, 'jp1')
    t = sub_once(r'^(GL_EXP = dTdrj.*)$',
                 r'\1' '\nj = next(k for k in range(1, npointsm1) if dTdr_v[k] <= dTdrj < dTdr_v[k+1])\njp1 = j + 1',
                 t, 'GL_EXP')
    # nucleation-rate figure: the fixed axes of the original hide the elasto-plastic curve -> autoscale
    m = re.search(r"^ax\d+\.set_ylabel\(r'Nucleation Sites", t, flags=re.M)
    if m:
        # elastic run of elastic_rate(): writes the rate and stops before the rate figure
        f0 = list(re.finditer(r'^fig\d+\s*,\s*ax\d+\s*=\s*plt\.subplots', t[:m.start()], flags=re.M))[-1].start()
        dump = ("if __import__('os').environ.get('GM_DUMP_RATE'):   # elastic rate for the comparison (figure_layout.elastic_rate)\n"
                "    __import__('numpy').savez(__import__('os').environ['GM_DUMP_RATE'], DT=DT_v[1:npointsm1], "
                "G=dTdr_v[1:npointsm1], Ihom=stat_hom_2nd_v[1:npointsm1], Ihet=stat_het_2nd_v[1:npointsm1], "
                "GT=GT_het_2nd_v[1:npointsm1], gam=gam_hom_v[1:npointsm1], tau=sigma_v[1:npointsm1])\n"
                "    raise SystemExit\n")
        t = t[:f0] + dump + t[f0:]
        m = re.search(r"^ax\d+\.set_ylabel\(r'Nucleation Sites", t, flags=re.M)
        k = t.find('\nplt.axis(', m.end())
        nxt = re.search(r'^fig\d+\s*,\s*ax\d+\s*=\s*plt\.subplots', t[m.end():], flags=re.M)
        if k >= 0 and (nxt is None or k < m.end() + nxt.start()):
            e = t.find('\n', k + 1)
            old = t[k + 1:e]
            # the texts of the original are in data coordinates of the original limits -> inside the axes
            ax = re.findall(r'^(ax\d+)\.set_ylabel\(', t[:m.end()], flags=re.M)[-1]
            seg = t[m.end():k]
            new_seg = seg
            t = t[:m.end()] + new_seg + t[k:]
            k += len(new_seg) - len(seg)
            e = t.find('\n', k + 1)
            overlay = (
                '# elastic rate of the same script (separate process) and end of the elastic range, for comparison\n'
                "_ov = elastic_rate(__file__) if (gm_surf.model == 'elastoplastic' and not __import__('os').environ.get('GM_DUMP_RATE')) else None\n"
                f'_n_ep = len({ax}.get_lines())\n'
                'if _ov is not None:\n'
                f"    {ax}.plot(_ov['DT'], _ov['Ihom'], color='#2a78d6', linewidth=1.5)\n"
                f"    {ax}.plot(_ov['DT'], _ov['Ihet'], color='#2a78d6', linewidth=2.0, linestyle='-.')\n"
                '_DTy = None\n'
                'if gm_surf.model == \'elastoplastic\' and gm_surf.yielded_at:\n'
                "    _iy = int(__import__('numpy').argmin(__import__('numpy').abs(rc_v[1:npointsm1] - gm_surf.yielded_at)))\n"
                '    _DTy = DT_v[1:npointsm1][_iy]\n'
                f"    {ax}.axvspan(DT_v[1], _DTy, color='#2a78d6', alpha=0.08, lw=0)\n"
                f"    {ax}.axvline(_DTy, color='k', linewidth=1.0, linestyle='-.')\n"
                'if gm_surf.model == \'elastoplastic\':\n'
                f'    _old = [_t.get_text() for _t in {ax}.get_legend().get_texts()] if {ax}.get_legend() else []\n'
                "    _lab = [('Homogeneous' if _i == 0 else 'Heterogeneous') + ', elasto-plastic' for _i in range(min(_n_ep, len(_old)))]\n"
                '    if _ov is not None:\n'
                "        _lab += ['Homogeneous, elastic', 'Heterogeneous, elastic']\n"
                '    if _DTy is not None:\n'
                "        _lab += [r'Yield (elastic range shaded), $\\Delta T$ = %.3f K' % _DTy]\n"
                f"    {ax}.legend({ax}.get_lines()[:len(_lab)], _lab, title='2nd-order statistics', frameon=True, edgecolor='k', shadow=True)\n")
            # the original may hold alloy names of other scripts (one of them outside its axes): keep this alloy's
            key = ALLOY_TEXT.get(alloy)
            alloy_text = (f"_al = [_tx for _tx in {ax}.texts if 'Alloy' in _tx.get_text()]   # one alloy name only\n"
                          'if len(_al) > 1:\n'
                          f"    _keep = [_tx for _tx in _al if {key!r} and {key!r} in _tx.get_text()] or _al[-1:]\n"
                          '    for _tx in _al:\n'
                          '        if _tx not in _keep:\n'
                          '            _tx.remove()\n')
            t = t[:k + 1] + ('#' + old + '   # original limits (they hide the elasto-plastic rate)\n'
                             'plt.autoscale(enable=True, axis="both")\n' + overlay + alloy_text +
                             f'for _k, _txt in enumerate(sorted({ax}.texts, key=lambda tx: -tx.get_position()[1])):   # texts inside the axes\n'
                             f'    _txt.set_transform({ax}.transAxes)\n'
                             '    _txt.set_position((0.27, 0.94 - 0.062 * _k))\n'
                             '    _txt.set_fontsize(min(_txt.get_fontsize(), 11))\n'
                             f'_leg = {ax}.get_legend()   # legend inside the axes, below the texts\n'
                             'if _leg is not None:\n'
                             "    _leg.set_loc('upper right')\n"
                             f'    _leg.set_bbox_to_anchor((1.0, 0.95 - 0.062 * len({ax}.texts)), transform={ax}.transAxes)\n'
                             '    for _lt in _leg.get_texts():\n'
                             '        _lt.set_fontsize(9)') + t[e:]
    # new figure: nucleation rate against the thermal gradient (the figure against the undercooling is kept)
    # texts of the original outside the axes (data coordinates of the original sweep) -> inside the axes
    t, k = re.subn(r'^plt\.show\(\)',
                   lambda _m: RATE_VS_GRADIENT + 'texts_inside_all()   # figure_layout.py: texts inside the axes\nplt.show()',
                   t, flags=re.M)
    if k:
        t = sub_once(r'^(from gm_elastoplastic import .*)$', r'\1\nfrom figure_layout import texts_inside_all, elastic_rate, vs_gradient_figure', t, 'figure_layout')
    header = ('"""\nELASTO-PLASTIC COPY (modelos_elasto-plastos). Original: Sim Luane Final/%s\n'
              'The Gurtin-Murdoch surface of the main loop is gm_elastoplastic.GMElastoPlasticSurface;\n'
              'GM_SURFACE_MODEL=elastic reproduces the original. See ../README.md.\n"""\n')
    return header, t


def main():
    os.makedirs(DST, exist_ok=True)
    for h in HELPERS:
        shutil.copy2(os.path.join(SRC, h), DST)
    for rel, alloy in MODELS.items():
        src = os.path.join(SRC, rel)
        text = open(src, encoding='utf-8').read()
        header, new = convert(text, alloy)
        base = os.path.basename(rel)[:-3].replace(' ', '_') + '_EP.py'
        with open(os.path.join(DST, base), 'w', encoding='utf-8') as fh:
            if new.startswith('#!'):
                first, rest = new.split('\n', 1)
                fh.write(first + '\n' + header % rel.replace('\\', '/') + rest)
            else:
                fh.write(header % rel.replace('\\', '/') + new)
        print(f'converted: {rel}  ->  nucleacao/{base}  ({alloy}, sigma_y = {SIGMA_Y[alloy]:g} N/m)')
    nc = os.path.join(DST, 'nao_convertidos')
    os.makedirs(nc, exist_ok=True)
    for rel in NOT_CONVERTED:
        shutil.copy2(os.path.join(SRC, rel), nc)
        print(f'copied without conversion: {rel}  ->  nucleacao/nao_convertidos/')


if __name__ == '__main__':
    main()
