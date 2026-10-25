"""Builds the elasto-plastic surface paper (Word, native equations)."""
import os
import re
import copy
import numpy as np
import latex2mathml.converter as l2m
from lxml import etree
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_TAB_ALIGNMENT
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = r'C:\Users\ileao\OneDrive\Documentos\Paper Nucleation\Manuscript_ElastoPlastic_Surface_NonEquilibrium_Nucleation_v7_final.docx'
XSLT = etree.XSLT(etree.parse(r'C:\Program Files\Microsoft Office\root\Office16\MML2OMML.XSL'))

doc = Document()
sec = doc.sections[0]
sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
for side in ('left_margin', 'right_margin'):
    setattr(sec, side, Cm(2.5))
sec.top_margin = sec.bottom_margin = Cm(2.5)
TEXT_W = 16.0  # cm

st = doc.styles['Normal']
st.font.name = 'Times New Roman'
st.font.size = Pt(12)
st.element.rPr.rFonts.set(qn('w:eastAsia'), 'Times New Roman')
st.paragraph_format.space_after = Pt(6)
st.paragraph_format.line_spacing = 1.5
for name, size in (('Heading 1', 13), ('Heading 2', 12), ('Title', 16)):
    h = doc.styles[name]
    h.font.name = 'Times New Roman'
    h.font.size = Pt(size)
    h.font.bold = True
    h.font.italic = False
    h.font.color.rgb = RGBColor(0, 0, 0)
    rpr = h.element.get_or_add_rPr()
    rpr.get_or_add_rFonts().set(qn('w:ascii'), 'Times New Roman')
    rpr.get_or_add_rFonts().set(qn('w:hAnsi'), 'Times New Roman')
    h.paragraph_format.space_before = Pt(12)
    h.paragraph_format.space_after = Pt(6)
    h.paragraph_format.line_spacing = 1.15

eq_no = [0]


MNS = 'http://schemas.openxmlformats.org/officeDocument/2006/math'


def omml(latex):
    # spaces of the converter are dropped by MML2OMML: use non-breaking spaces
    latex = latex.replace(r'\qquad', r'\text{\u00a0\u00a0\u00a0\u00a0}').replace(r'\quad', r'\text{\u00a0\u00a0}')
    latex = latex.replace(r'\text{\u00a0\u00a0\u00a0\u00a0}', '\\text{\u00a0\u00a0\u00a0\u00a0}').replace(r'\text{\u00a0\u00a0}', '\\text{\u00a0\u00a0}')
    mml = l2m.convert(latex)
    root = XSLT(etree.fromstring(mml)).getroot()
    # \dot -> limUpp with U+02D9; Word expects an accent (m:acc, U+0307)
    for lu in root.iter(f'{{{MNS}}}limUpp'):
        lim = lu.find(f'{{{MNS}}}lim')
        if lim is not None and ''.join(lim.itertext()) in ('\u02d9', '\u00a8'):
            acc = etree.Element(f'{{{MNS}}}acc')
            pr = etree.SubElement(acc, f'{{{MNS}}}accPr')
            ch = etree.SubElement(pr, f'{{{MNS}}}chr')
            ch.set(f'{{{MNS}}}val', '\u0307' if ''.join(lim.itertext()) == '\u02d9' else '\u0308')
            acc.append(lu.find(f'{{{MNS}}}e'))
            lu.getparent().replace(lu, acc)
    return root


def add_math(par, latex):
    el = omml(latex)
    par._p.append(etree.fromstring(etree.tostring(el)) if False else _to_oxml(el))


def _to_oxml(el):
    from docx.oxml import parse_xml
    return parse_xml(etree.tostring(el))


def rich(par, text, bold=False, italic=False, size=None):
    """Text with inline LaTeX between $...$."""
    for i, part in enumerate(re.split(r'\$([^$]+)\$', text)):
        if i % 2:
            par._p.append(_to_oxml(omml(part)))
        elif part:
            r = par.add_run(part)
            r.bold, r.italic = bold, italic
            if size:
                r.font.size = Pt(size)
    return par


def P(text, align=WD_ALIGN_PARAGRAPH.JUSTIFY, **kw):
    par = doc.add_paragraph()
    par.alignment = align
    return rich(par, text, **kw)


def H1(t):
    doc.add_heading(t, level=1)


def H2(t):
    doc.add_heading(t, level=2)


def EQ(latex, label=None):
    eq_no[0] += 1
    par = doc.add_paragraph()
    pf = par.paragraph_format
    pf.tab_stops.add_tab_stop(Cm(TEXT_W / 2), WD_TAB_ALIGNMENT.CENTER)
    pf.tab_stops.add_tab_stop(Cm(TEXT_W), WD_TAB_ALIGNMENT.RIGHT)
    pf.space_before = Pt(3)
    pf.space_after = Pt(3)
    par.add_run('\t')
    par._p.append(_to_oxml(omml(latex)))
    par.add_run(f'\t({eq_no[0]})')
    if label:
        EQS[label] = eq_no[0]
    return eq_no[0]


EQS = {}


def FIG(path, caption, width=TEXT_W):
    par = doc.add_paragraph()
    par.alignment = WD_ALIGN_PARAGRAPH.CENTER
    par.add_run().add_picture(path, width=Cm(width))
    cap = doc.add_paragraph()
    cap.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    cap.paragraph_format.line_spacing = 1.15
    caption_runs(cap, caption)


def caption_runs(cap, caption):
    label, rest = caption.split('. ', 1)
    r = cap.add_run(label + '. ')
    r.bold = True
    r.font.size = Pt(10)
    rich(cap, rest, size=10)


def shade(cell, hexcolor):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement('w:shd')
    shd.set(qn('w:val'), 'clear')
    shd.set(qn('w:color'), 'auto')
    shd.set(qn('w:fill'), hexcolor)
    tcPr.append(shd)


def TABLE(caption, header, rows, widths):
    cap = doc.add_paragraph()
    cap.paragraph_format.line_spacing = 1.15
    cap.paragraph_format.keep_with_next = True
    caption_runs(cap, caption)
    t = doc.add_table(rows=1 + len(rows), cols=len(header))
    t.style = 'Table Grid'
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for j, w in enumerate(widths):
        for row in t.rows:
            row.cells[j].width = Cm(w)
    for j, h in enumerate(header):
        c = t.rows[0].cells[j]
        c.paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.CENTER
        c.paragraphs[0].paragraph_format.line_spacing = 1.0
        rich(c.paragraphs[0], h, bold=True, size=9)
        shade(c, 'E7E6E6')
    for i, r in enumerate(rows):
        for j, v in enumerate(r):
            c = t.rows[i + 1].cells[j]
            c.paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.LEFT if j == 0 else WD_ALIGN_PARAGRAPH.CENTER
            c.paragraphs[0].paragraph_format.line_spacing = 1.0
            c.paragraphs[0].paragraph_format.space_after = Pt(0)
            rich(c.paragraphs[0], v, size=9)
    doc.add_paragraph()


# ------------------------------------------------------------------ title
t = doc.add_paragraph()
t.alignment = WD_ALIGN_PARAGRAPH.CENTER
t.paragraph_format.line_spacing = 1.15
r = t.add_run('Elasto-Plastic Gurtin–Murdoch Surface in Non-Equilibrium Phase Nucleation: '
              'Stored Permanent Surface Deformation, Elastic Restoration and the Thermal-Field Tensor')
r.bold = True
r.font.size = Pt(15)
P('I. L. Ferreira(*)', align=WD_ALIGN_PARAGRAPH.CENTER)
P('Faculty of Mechanical Engineering, Federal University of Pará, UFPA, Augusto Corrêa Avenue 1, 66075-110, '
  'Belém, PA, Brazil', align=WD_ALIGN_PARAGRAPH.CENTER, size=11)
P('*Corresponding author email: ileao@ufpa.br', align=WD_ALIGN_PARAGRAPH.CENTER, size=11)

H1('Abstract')
P('In previous works, the Gibbs–Thomson coefficient was generalized as a thermal-field tensor derived from the first '
  'and second laws of thermodynamics. Its second-order solution, based on Gurtin–Murdoch elastic surface-stress '
  'theory, predicted unbounded excess surface stress at increasing thermal gradients. Here, the model is extended to '
  'elasto-plasticity using non-equilibrium thermodynamics with internal variables. The areal strain is decomposed '
  'into elastic and permanent components, and the surface Helmholtz free energy, Clausius–Duhem inequality, areal '
  'yield condition, and maximum dissipation principle determine permanent strain evolution, plastic work, and '
  'dissipation. Permanent deformation stores energy available for phase nucleation, while plastic-flow dissipation '
  'releases heat. A linear relaxation law describes elastic restoration driven by phase-transformation density '
  'changes. These corrections are incorporated into first- and second-order thermal-field tensors and the nucleation '
  'rate of the FCC_A1 phase in Al–3Cu–5Nb–0.1Fe, Al–0.8Si–0.6Mg–0.2Fe, and Al–3Cu–0.5Mg–0.2Fe alloys. Interfacial '
  'yielding occurs at gradients of 560–925 K m⁻¹. Thereafter, the elasto-plastic tensor remains 0.6–4.2% below its '
  'elastic counterpart, while nucleation-rate maxima decrease by 4–19%. Surface tension differs from the elastic '
  'prediction by 0.6–5.3%. Across gradients up to 10⁶ K m⁻¹, surface energy reaches 19–25 times its equilibrium '
  'value at 10⁵ K m⁻¹.')
P('Keywords: Phase nucleation; Non-equilibrium thermodynamics; Surface plasticity; Gurtin–Murdoch surface stress; '
  'Thermal-field tensor; Gibbs–Thomson coefficient; Aluminum alloys.', size=11)

# ------------------------------------------------------------------ 1
H1('1. Introduction')
P('Nucleation during solidification takes place far from equilibrium: the experimental cooling rates of transient '
  'solidification are one to three orders of magnitude higher than those at which equilibrium is reached, and the '
  'undercooling and the solid–liquid interfacial energy depend on the local thermal gradient [1–4]. In a series of works, '
  'the author proposed a nucleation formulation written in terms of primitive thermodynamic variables, in which the first '
  'law defines a thermal field whose gradient, $\\nabla T$, displaces the solid–liquid equilibrium, and the Gibbs–Thomson '
  'coefficient becomes a thermal-field tensor, $\\Gamma=A\\cdot\\nabla T$, where $A$ is the newly created or deformed surface area '
  'normal to the thermal field [1, 2]. The maximization of the total Gibbs free energy with respect to the nucleus radius '
  'gives a quadratic equation, whose full second-order solution provides the critical radius and the tensor, while the '
  'first-order solution is an approximation of it [1, 2, 5].')
P('The solution requires the surface stress, surface tension and surface energy of the nucleus, which were computed '
  'with the surface elasticity of Gurtin and Murdoch [6, 7] in the continuum-mechanics domain. In that description the '
  'interface is a two-dimensional elastic body with a residual (compressive) tension $\\sigma_0$ and surface Lamé moduli '
  '$\\lambda_0$ and $\\mu_0$, and the surface energy grows as the radius decreases from the equilibrium reference '
  '$r_{eq}$. The previous calculations assumed explicitly that the nucleation generates a homogeneous surface stress and '
  'that no plastic deformation occurs [2]. Under this assumption the elastic excess stress grows without bound when the '
  'thermal gradient increases and the critical radius decreases, so that the surface energy reaches values several times '
  'larger than the equilibrium value $\\gamma_0$ at gradients of $10^4$–$10^5$ K m⁻¹ [2, 3].')
P('A real interface cannot sustain an unlimited areal stress: beyond a threshold, the excess of atoms or sites created '
  'by the areal change is accommodated permanently, by the same mechanisms that produce plastic flow in the bulk. In bulk '
  'metals, the work of plastic deformation is partly dissipated as heat and partly stored as defect energy [8], and this '
  'stored energy is the driving force for recovery, recrystallization and phase transformations. The thermodynamic '
  'description of such irreversible processes requires internal variables and the second law in the form of the '
  'Clausius–Duhem inequality [9–12]. In the thermal-field formulation, moreover, the path-dependent work $W$ was '
  'introduced from the first law but neglected in the elastic solution [2]; the plastic work of the interface is the '
  'natural content of that term.')
P('The objective of the present work is to generalize the Gurtin–Murdoch surface to elasto-plasticity within '
  'non-equilibrium thermodynamics, and to introduce all resulting corrections of the surface energy into the '
  'thermal-field tensor and into the nucleation rate. The corrections are: (i) an areal yield condition that limits the '
  'excess surface stress and marks the end of the elastic range; (ii) the permanent areal strain; (iii) the energy of '
  'the permanent surface deformation, stored in the interface and available for nucleation; (iv) the hardening energy '
  'and the heat dissipated by the plastic flow; (v) the consistent elasto-plastic tangent in the radial derivative of the surface energy, which enters '
  'both solutions of the quadratic equation; and (vi) the elastic restoration of the recoverable part of the surface '
  'energy, driven by the rate of density change of the transformation, which leaves the stored energy of the permanent '
  'deformation in the interface. The analysis of the secondary dendrite arm spacing (SDAS) with the '
  'elasto-plastic tensor is outside the scope of this paper and will be reported separately.')

# ------------------------------------------------------------------ 2
H1('2. Theoretical Background')
H2('2.1. Thermal field and thermal-field tensor')
P('Following Refs. [1, 2], the first law, $Q-W=E$, defines a thermal field gradient that can be split into a '
  'path-independent part, associated with the energy $E$, and a path-dependent part, associated with the work $W$:')
EQ(r'\nabla T\cong\nabla\left(\frac{\partial T}{\partial E}\right)E+\nabla\left(\frac{\partial T}{\delta W}\right)W', 'grad')
P('and the internal energy of the interface contains the surface term $A_0\\,\\Sigma:\\varepsilon-\\gamma\\,\\delta:\\varepsilon$ [2]. The '
  'total Gibbs free energy of a nucleus of radius $r_C$ and contact angle $\\theta$ is')
EQ(r'\Delta G=\frac{1}{3}\pi r_C^3\,\Delta S_V\Delta T+\pi r_C^2\,\gamma_{SL}f(\theta)=\frac{1}{3}\pi r_C^3\,\Delta G_V+\pi r_C^2\,\Delta G_S f(\theta)', 'dG')
P('with $f(\\theta)=2-3\\cos\\theta+\\cos^3\\theta$ and $\\Delta G_S=\\gamma_{SL}$. Since $\\Delta S_V$, $\\Delta T$, $\\gamma_{SL}$ and $f(\\theta)$ '
  'depend on the radius, $\\partial\\Delta G/\\partial r=0$ gives a quadratic equation whose full second-order solution and '
  'first-order approximation are [2, 5]')
EQ(r'\Gamma^{2nd}=-\frac{3}{4}\,\frac{\Delta G_V+\frac{\partial\Delta G_S}{\partial r}+\Delta G_C}{\frac{\partial\Delta S_V}{\partial r}}=A\cdot\nabla T,\qquad \Gamma^{1st}=-\frac{\Delta G_S}{\Delta S_V+\frac{\partial\Delta S_S}{\partial r}+\Delta S_C}', 'GT')
P('where $\\Delta G_C$ is the configurational free energy and $\\Delta S_S=\\Delta T^{-1}\\partial\\gamma_{SL}/\\partial r$. For the '
  'homogeneous nucleus ($\\theta=\\pi$, $\\Delta G_C=0$) the undercooling and the thermal-field gradient at the critical '
  'radius are')
EQ(r'\Delta T=\frac{2\gamma_{SL}}{r\left|\Delta S_V\right|},\qquad \nabla T=\frac{\Delta T}{8\pi r},\qquad \Gamma^{2nd}=4\pi r^2\,\nabla T=\frac{r\Delta T}{2}', 'DT')
P('and the first-order tensor can be written as a correction of the second-order one by the logarithmic derivative of '
  'the surface energy:')
EQ(r'\Gamma^{1st}=-\frac{\gamma_{SL}}{\Delta S_V+\frac{1}{\Delta T}\frac{\partial\gamma_{SL}}{\partial r}}=\frac{\Gamma^{2nd}}{1+\frac{r}{2\gamma_{SL}}\left|\frac{\partial\gamma_{SL}}{\partial r}\right|}', 'G1')
P(f'Equation ({EQS["G1"]}) shows that every correction of $\\gamma_{{SL}}$ and of its radial derivative enters both '
  'solutions; the second-order solution is the complete one and the first-order solution contains the derivative of the '
  'surface energy as an implicit reduction factor.')

H2('2.2. Surface free energy from the metric of the sphere: variation with the area')
P('The interface is a two-dimensional surface with induced metric $ds^2=g_{ij}\\,dx^i dx^j$. For the surface of a '
  'nucleus of radius $r$, in the coordinates $(\\theta,\\phi)$, $g_{\\theta\\theta}=r^2$, $g_{\\phi\\phi}=r^2\\sin^2\\theta$, and the '
  'invariant area element and the area follow from the metric alone:')
EQ(r'dA=\sqrt{\det\left(g_{ij}\right)}\,d\theta\,d\phi=r^2\sin\theta\,d\theta\,d\phi,\qquad A=\int_0^{2\pi}\int_0^{\pi}r^2\sin\theta\,d\theta\,d\phi=4\pi r^2,\qquad V(A)=\frac{A^{3/2}}{6\sqrt{\pi}}', 'area')
P('The fundamental geometric variable of the interface is therefore the area $A$, defined by the metric, and not the '
  'radius, which only parametrizes $A$ in the spherical case. The surface free energy is the functional')
EQ(r'G_S=\int_S\gamma\,dA,\qquad \frac{dG_S}{dA}=\gamma+A\frac{d\gamma}{dA}', 'GSA')
P('so that, when the surface energy depends on the area or on the curvature, the variation of $G_S$ is not simply '
  '$\\gamma$. For a uniform $\\gamma$ and a constant $\\Delta g_V$, $\\Delta G(A)=\\Delta g_V V(A)+\\gamma A$, and the critical '
  'condition taken with respect to the area,')
EQ(r'\frac{d\Delta G}{dA}=\frac{\Delta g_V}{4\sqrt{\pi}}A^{1/2}+\gamma=\frac{\Delta g_V\,r}{2}+\gamma=0\quad\Rightarrow\quad r_C=-\frac{2\gamma}{\Delta g_V}', 'CNTA')
P('recovers the classical critical radius (CNT). For a sphere that remains spherical, varying $A$ or $r$ is equivalent '
  'by the chain rule, $dA=8\\pi r\\,dr$; the change of variable alone does not modify the CNT. The modification comes from '
  'the dependence of $\\gamma$ on the area, through the metric, the curvature and the surface stress: in Eq. '
  f'({EQS["G1"]}), $\\frac{{r}}{{2\\gamma_{{SL}}}}\\left|\\frac{{\\partial\\gamma_{{SL}}}}{{\\partial r}}\\right|=\\frac{{A}}{{\\gamma_{{SL}}}}\\left|\\frac{{\\partial\\gamma_{{SL}}}}{{\\partial A}}\\right|$, '
  'i.e., the first-order tensor is reduced by the logarithmic derivative of the surface energy with respect to the '
  'area. The same variable governs the deformation of the interface: the areal strain used below, '
  '$E_A=2\\ln\\left[y(r)/y(r_{eq})\\right]$, tends to $\\ln(A/A_{eq})$ when $(3\\lambda+2\\mu)r\\gg 2(\\lambda_0+\\mu_0)$, and the '
  'elasto-plastic flow of Section 3 is a flow of area.')

H2('2.3. Elastic Gurtin–Murdoch surface and surface free energy')
P('The solution of the thermal-field tensor requires the solution of the surface-stress tensor $\\Sigma$ of Gurtin and '
  'Murdoch [6, 7], which satisfies the surface equilibrium and the constitutive relation')
EQ(r'\nabla_S\cdot\Sigma=\Theta\,n+\rho_0\ddot u,\qquad \Sigma=\sigma_0 I+2\left(\mu_0-\sigma_0\right)E+\left(\lambda_0-\sigma_0\right)\left(\mathrm{tr}\,E\right)I+\sigma_0\nabla_S u,\qquad E=\frac{1}{2}\left(Du+Du^T\right)', 'GM')
P('where $\\Theta$, $E$ and $u$ are the stress, strain and displacement fields. On the nucleus surface, in spherical '
  'coordinates, the metric tensor and the contravariant metric of the initial (state 0) configuration are')
EQ(r'g_{ij}=\left(\begin{array}{ccc}1&0&0\\0&r^2&0\\0&0&r^2\sin^2\theta\end{array}\right),\qquad g_0^{ij}=\left(\begin{array}{ccc}1&0&0\\0&r_0^{-2}&0\\0&0&r_0^{-2}\sin^{-2}\theta_0\end{array}\right)', 'metric')
P('Four formulations are written in terms of the metric, of the Ricci scalar $R$ and of the work of the surface stress '
  'on the area [2]. The surface free energy:')
EQ(r'\Delta G_S=\gamma_0\,\Phi\left(r_0^2R\right)-\frac{1}{2A}\int_S\Sigma^{ij}g_{ij}\,dA', 'GSgen')
P('the surface free-energy potential:')
EQ(r'\frac{\partial\Delta G_S}{\partial r}=\partial_k\left[\gamma_0\,\Phi\left(r_0^2R\right)-\frac{1}{2A}\int_S\Sigma^{ij}g_{ij}\,dA\right]=\Delta G_{S0}\left(1+r_0^2R\right)\frac{g_0^{ij}\nabla_k g_{ij}}{\left(g_{mn}g_0^{mn}\right)^2}-\frac{1}{2A}\partial_k\int_S\Sigma^{ij}g_{ij}\,dA', 'dGSgen')
P('the surface entropy:')
EQ(r'\Delta S_S=\frac{\gamma_0}{\Delta T}\Phi\left(r_0^2R\right)-\frac{1}{2A\,\Delta T}\int_S\Sigma^{ij}g_{ij}\,dA', 'SSgen')
P('and the surface entropy potential:')
EQ(r'\frac{\partial\Delta S_S}{\partial r}=\frac{1}{\Delta T}\partial_k\left[\gamma_0\,\Phi\left(r_0^2R\right)-\frac{1}{2A}\int_S\Sigma^{ij}g_{ij}\,dA\right]=-\frac{\Delta G_{S0}}{\Delta T}\left(1+r_0^2R\right)\frac{g_0^{ij}\nabla_k g_{ij}}{\left(g_{mn}g_0^{mn}\right)^2}-\frac{1}{\Delta T\,2A}\partial_k\int_S\Sigma^{ij}g_{ij}\,dA', 'dSSgen')
P('where $\\Phi$ is a dimensionless geometric function whose argument, $r_0^2R$, is the only admissible scalar curvature '
  'invariant, $\\Delta G_{S0}=\\gamma_0$ is the surface energy of state 0 and $r_0$ its radius (the equilibrium radius, '
  '$r_0=r_{eq}$). For a 2-sphere, $R=2/r_0^2$, $\\nabla_k R=0$, $\\nabla_k g_0^{ij}=0$ and $k=r$, and the surface free '
  'energy and its potential, Eqs. '
  f'({EQS["GSgen"]}) and ({EQS["dGSgen"]}), become')
EQ(r'\Delta G_S^{sph}=\Delta G_{S0}\frac{3}{1+\frac{r^2}{r_0^2}+\frac{r^2\sin^2\theta}{r_0^2\sin^2\theta_0}}-\frac{1}{4\pi r_0^2}\int_{A_0}^{A_1}\Sigma\,dA', 'GSsph')
EQ(r'\left.\frac{\partial\Delta G_S}{\partial r}\right|^{sph}=\Delta G_{S0}\frac{-6\left[\frac{r\sin^2\theta}{r_0^2\sin^2\theta_0}+\frac{r}{r_0^2}\right]}{\left(1+\frac{r^2}{r_0^2}+\frac{r^2\sin^2\theta}{r_0^2\sin^2\theta_0}\right)^2}-\frac{1}{\frac{r_0^2}{r^2}}\Sigma_r(r)', 'dGSsph')
P('In the elastic calculations of Refs. [1–3], the spherical solution of the bulk–surface equilibrium with bulk Lamé moduli '
  '$\\lambda$ and $\\mu$ reduces the work of the surface stress to the excess surface stress $s$ and the geometric function to '
  '$\\Phi=(r_{eq}/r)^2$, which is the particular form used in the numerical results of the present paper.')
P('Here $\\lambda_0$ and $\\mu_0$ are the surface Lamé moduli and $\\sigma_0$ the residual surface tension, which is '
  'compressive ($\\sigma_0<0$). For the spherical nucleus, the solution '
  'of the bulk–surface equilibrium with bulk Lamé moduli $\\lambda$ and $\\mu$ gives the areal (logarithmic) strain of the '
  'surface relative to the equilibrium reference $r_{eq}$ and the excess surface stress [1, 2]:')
EQ(r'E_A=2\ln\frac{y(r)}{y(r_{eq})},\qquad y(r)=\left(3\lambda+2\mu\right)r+2\left(\lambda_0+\mu_0\right),\qquad s=C_A E_A,\quad C_A=-\sigma_0=\left|\sigma_0\right|', 'EA')
P('so that $E_A<0$ and $s<0$ for $r<r_{eq}$. The surface tension and the surface energy are')
EQ(r'\tau=\sigma_0+s,\qquad \gamma_{SL}^{el}=\gamma_0\left(\frac{r_{eq}}{r}\right)^2-s=\gamma_0\left(\frac{r_{eq}}{r}\right)^2+\sigma_0 E_A', 'gel')
P('where $\\gamma_0$ is the equilibrium surface energy. Both terms of Eq. '
  f'({EQS["gel"]}) are positive: the first is the geometric (curvature) term and the second, $\\left|\\sigma_0\\right|\\left|E_A\\right|$, '
  'is the work of the residual tension on the change of area. The radial derivative used in Eqs. '
  f'({EQS["GT"]})–({EQS["G1"]}) is')
EQ(r'\frac{\partial\gamma_{SL}^{el}}{\partial r}=-\frac{2\gamma_0 r_{eq}^2}{r^3}+\sigma_0\frac{\partial E_A}{\partial r},\qquad \frac{\partial E_A}{\partial r}=\frac{2}{r\left(1+\alpha_r\right)},\quad \alpha_r=\frac{2\left(\left|\sigma_0\right|+2\left(\lambda_0+\mu_0\right)\right)}{r\left(3\lambda+2\mu\right)}', 'dgel')

# ------------------------------------------------------------------ 3
H1('3. Elasto-Plastic Surface in Non-Equilibrium Thermodynamics')
H2('3.1. Kinematics and internal variables')
P('The areal strain is split additively into an elastic (recoverable) part and a permanent part, and an accumulated '
  'areal plastic strain $\\alpha$ is introduced as an internal variable of the interface [9–12]:')
EQ(r'E_A=E_A^e+E_A^p,\qquad \dot{\alpha}=\left|\dot E_A^p\right|,\qquad \alpha\ge 0', 'split')
P('The permanent areal strain is the part of the areal change that is not recovered on unloading. It could be released '
  'by moving the equilibrium reference radius of the nucleus,')
EQ(r'A_{eq,p}=A_{eq}\exp\left(E_A^p\right)\quad\Rightarrow\quad r_{eq,p}=r_{eq}\exp\left(\frac{E_A^p}{2}\right)', 'ref')
P('but in the present formulation it is not released: the permanent area remains stored in the interface and the '
  'equilibrium reference remains $r_{eq}$ (Section 3.5).')

H2('3.2. Surface free energy and Clausius–Duhem inequality')
P('The surface Helmholtz free energy per unit reference area is taken as the sum of a recoverable (elastic) part and '
  'a hardening (stored defect) part:')
EQ(r'\psi^s\left(E_A^e,\alpha\right)=\psi_0^s+\frac{1}{2}C_A\left(E_A^e\right)^2+\psi_p(\alpha),\qquad \psi_p=\frac{1}{2}H\alpha^2', 'psi')
P('where $H$ is the areal hardening modulus. For the isothermal interface, the second law in the form of the '
  'Clausius–Duhem inequality of the surface is [9, 10, 13]')
EQ(r'\mathcal{D}^s=s\,\dot E_A-\dot\psi^s=\left(s-\frac{\partial\psi^s}{\partial E_A^e}\right)\dot E_A^e+s\,\dot E_A^p-q\,\dot\alpha\ge 0,\qquad q=\frac{\partial\psi^s}{\partial\alpha}=H\alpha', 'CD')
P('Since the inequality must hold for arbitrary elastic processes, the Coleman–Noll procedure [9] gives the state law '
  '$s=\\partial\\psi^s/\\partial E_A^e=C_A E_A^e$ and the reduced dissipation')
EQ(r'\mathcal{D}^s=s\,\dot E_A^p-q\,\dot\alpha\ge 0,\qquad \Pi^s=\frac{\mathcal{D}^s}{T}\ge 0', 'Dred')
P('where $\\Pi^s$ is the entropy production of the interface per unit area. The surface remains in a non-equilibrium '
  'state as long as $\\mathcal{D}^s>0$.')

H2('3.3. Areal yield condition, flow rule and return mapping')
P('The elastic domain of the interface is bounded by an areal yield condition on the excess surface stress, with an '
  'areal yield stress $\\sigma_y^s$ [N m⁻¹]:')
EQ(r'\phi\left(s,q\right)=\left|s\right|-\left(\sigma_y^s+q\right)\le 0', 'yield')
P('The yield point, $\\phi=0$, marks the end of the elastic range of the interface and the beginning of the plastic one. '
  'On a spherical nucleus the two principal surface stresses are equal, so that the deviatoric (von Mises) part of the '
  'surface stress vanishes and only the areal (isotropic) mode can yield. The principle of maximum dissipation [12, 14] '
  'gives the associative flow rule and the Kuhn–Tucker conditions')
EQ(r'\dot E_A^p=\dot\lambda\,\mathrm{sign}(s),\qquad \dot\alpha=\dot\lambda,\qquad \dot\lambda\ge 0,\quad \phi\le 0,\quad \dot\lambda\,\phi=0', 'flow')
P('with which the reduced dissipation becomes $\\mathcal{D}^s=\\sigma_y^s\\dot\\lambda\\ge 0$, so that the second law is '
  'satisfied identically. As the radius decreases from $r_{eq}$, Eqs. '
  f'({EQS["EA"]}) and ({EQS["yield"]})–({EQS["flow"]}) are integrated by the return mapping [12]:')
EQ(r's^{tr}=C_A\left(E_{A,n+1}-E_{A,n}^p\right),\quad \phi^{tr}=\left|s^{tr}\right|-\left(\sigma_y^s+H\alpha_n\right),\quad \Delta\lambda=\frac{\left\langle\phi^{tr}\right\rangle}{C_A+H}', 'rm1')
EQ(r's_{n+1}=s^{tr}-C_A\Delta\lambda\,\mathrm{sign}\left(s^{tr}\right),\qquad E_{A,n+1}^p=E_{A,n}^p+\Delta\lambda\,\mathrm{sign}\left(s^{tr}\right),\qquad \alpha_{n+1}=\alpha_n+\Delta\lambda', 'rm2')
P('where $\\langle x\\rangle=\\max(x,0)$. The consistent elasto-plastic tangent of the excess stress and the rate of the '
  'permanent strain during yielding are')
EQ(r'\frac{ds}{dE_A}=\frac{C_A H}{C_A+H},\qquad \frac{dE_A^p}{dE_A}=\frac{C_A}{C_A+H}', 'tangent')
P('For $H=0$ (perfect plasticity), the excess stress is limited to $s=-\\sigma_y^s$ and the elastic component of the surface tension to '
  '$\\tau=\\sigma_0-\\sigma_y^s$, independently of the radius, after the yield point.')

H2('3.4. Plastic work: stored energy and dissipation')
P('The plastic work of the interface per unit area is split into the stored hardening energy and the heat '
  'dissipated:')
EQ(r'w_p=\int s\,dE_A^p=\psi_p+q_{diss},\qquad \psi_p=\frac{1}{2}H\alpha^2,\qquad q_{diss}=\sigma_y^s\alpha', 'wp')
P(f'In terms of Eq. ({EQS["grad"]}), the dissipated part $q_{{diss}}$ is a path-dependent work and belongs to the term '
  '$\\nabla\\left(\\partial T/\\delta W\\right)W$ of the thermal field, which was neglected in the elastic solution, whereas the stored '
  'part enters the internal energy of the interface and therefore the path-independent term.')

H2('3.5. Surface energy of the elasto-plastic interface')
P('The surface energy that enters the nucleation formulation must contain both elastic and plastic terms. The '
  'hypothesis adopted here is that the permanent surface deformation complements the elastic one: the areal change of '
  'the nucleus is the same, $E_A=E_A^e+E_A^p$, but beyond the yield point it is accommodated permanently, and the energy '
  'associated with the permanent part is not lost but stored in the interface, available for the nucleation of '
  'phases once the thermodynamic conditions are favourable. In the same way as the elastic term of Eq. '
  f'({EQS["gel"]}), which is the work of the residual tension on the elastic areal change, the permanent areal change is '
  'also worked by the residual tension; only the heat dissipated by the plastic flow, Eq. '
  f'({EQS["wp"]}), leaves the interface:')
EQ(r'\gamma_{SL}^{ep}=\gamma_0\left(\frac{r_{eq}}{r}\right)^2+\left|\sigma_0\right|\left(\left|E_A^e\right|+\left|E_A^p\right|\right)-q_{diss}=\gamma_{SL}^{el}-\sigma_y^s\alpha', 'gep')
P('The terms are the geometric (curvature) term, the stored elastic energy, the stored energy of the permanent '
  'deformation and the dissipated heat. Since the stored energy of the permanent area equals the elastic energy that '
  'the same areal change would have, the elasto-plastic surface energy is the elastic one minus the heat produced by '
  'the plastic flow, i.e., minus the entropy produced in the interface multiplied by the temperature. It therefore grows '
  'with the thermal gradient as the elastic surface energy and remains below it only by $q_{diss}$. Two limiting '
  f'cases follow from Eq. ({EQS["gep"]}): without yielding ($\\sigma_y^s\\to\\infty$, $\\alpha=0$) and without yield '
  f'stress ($\\sigma_y^s\\to 0$, $q_{{diss}}\\to 0$) it reduces exactly to Eq. ({EQS["gel"]}); the plastic correction is '
  'largest at intermediate yield stresses.')
P('In the general notation of Eqs. '
  f'({EQS["GM"]})–({EQS["dGSsph"]}), the surface strain splits as $E=E^e+E^p$ and the surface stress depends on the elastic '
  'strain only, $\\Sigma=\\sigma_0 I+\\mathbb{L}:E^e+\\sigma_0\\nabla_S u$, with $\\mathbb{L}:E=2(\\mu_0-\\sigma_0)E+(\\lambda_0-\\sigma_0)'
  '(\\mathrm{tr}\\,E)I$. The permanent strain carries the stored stress $\\Sigma_p=\\mathbb{L}:E^p$, so that $\\Sigma+\\Sigma_p$ is the '
  'stress of the total (elastic and permanent) deformation, and the yield condition acts on the isotropic part of the '
  'excess stress, $s=\\frac{1}{2}g_{ij}\\left(\\Sigma^{ij}-\\sigma_0 g^{ij}\\right)$. The four formulations of Eqs. '
  f'({EQS["GSgen"]})–({EQS["dSSgen"]}) for the elasto-plastic interface are the surface free energy,')
EQ(r'\Delta G_S^{ep}=\gamma_0\,\Phi\left(r_0^2R\right)-\frac{1}{2A}\int_S\left(\Sigma^{ij}+\Sigma_p^{ij}\right)g_{ij}\,dA-\frac{1}{A}\int_S q_{diss}\,dA', 'GSep')
P('the surface free-energy potential,')
EQ(r'\frac{\partial\Delta G_S^{ep}}{\partial r}=\Delta G_{S0}\left(1+r_0^2R\right)\frac{g_0^{ij}\nabla_k g_{ij}}{\left(g_{mn}g_0^{mn}\right)^2}-\frac{1}{2A}\partial_k\int_S\left(\Sigma^{ij}+\Sigma_p^{ij}\right)g_{ij}\,dA-\sigma_y^s\,\partial_k\alpha', 'dGSep')
P('the surface entropy,')
EQ(r'\Delta S_S^{ep}=\frac{\gamma_0}{\Delta T}\Phi\left(r_0^2R\right)-\frac{1}{2A\,\Delta T}\int_S\left(\Sigma^{ij}+\Sigma_p^{ij}\right)g_{ij}\,dA-\frac{1}{A\,\Delta T}\int_S q_{diss}\,dA', 'SSep')
P('and the surface entropy potential,')
EQ(r'\frac{\partial\Delta S_S^{ep}}{\partial r}=-\frac{\Delta G_{S0}}{\Delta T}\left(1+r_0^2R\right)\frac{g_0^{ij}\nabla_k g_{ij}}{\left(g_{mn}g_0^{mn}\right)^2}-\frac{1}{\Delta T\,2A}\partial_k\int_S\left(\Sigma^{ij}+\Sigma_p^{ij}\right)g_{ij}\,dA-\frac{\sigma_y^s}{\Delta T}\partial_k\alpha', 'dSSep')
P('For the 2-sphere, the surface free energy and its potential, Eqs. '
  f'({EQS["GSep"]}) and ({EQS["dGSep"]}), become')
EQ(r'\Delta G_S^{ep,sph}=\Delta G_{S0}\frac{3}{1+\frac{r^2}{r_0^2}+\frac{r^2\sin^2\theta}{r_0^2\sin^2\theta_0}}-\frac{1}{4\pi r_0^2}\int_{A_0}^{A_1}\left(\Sigma+\Sigma_p\right)dA-\sigma_y^s\alpha', 'GSepsph')
EQ(r'\left.\frac{\partial\Delta G_S^{ep}}{\partial r}\right|^{sph}=\Delta G_{S0}\frac{-6\left[\frac{r\sin^2\theta}{r_0^2\sin^2\theta_0}+\frac{r}{r_0^2}\right]}{\left(1+\frac{r^2}{r_0^2}+\frac{r^2\sin^2\theta}{r_0^2\sin^2\theta_0}\right)^2}-\frac{1}{\frac{r_0^2}{r^2}}\left(\Sigma_r+\Sigma_{p,r}\right)(r)-\sigma_y^s\frac{\partial\alpha}{\partial r}', 'dGSepsph')
P('where, during yielding, $\\Sigma_r$ is computed with the consistent tangent of Eq. '
  f'({EQS["tangent"]}), $\\Sigma_{{p,r}}$ with $\\partial E_A^p/\\partial E_A=C_A/(C_A+H)$, and $\\partial\\alpha/\\partial r$ from the flow '
  f'rule. Without permanent strain ($E^p=0$, '
  f'$\\alpha=0$) Eqs. ({EQS["GSep"]})–({EQS["dGSepsph"]}) reduce to the elastic forms of Eqs. ({EQS["GSgen"]})–({EQS["dGSsph"]}). '
  f'Eq. ({EQS["gep"]}) is their particular form with $\\Phi=(r_{{eq}}/r)^2$ and the spherical work of the surface stress; with the '
  f'spherical function of Eq. ({EQS["GSepsph"]}) and $\\theta=\\theta_0$, the surface energy of the elasto-plastic interface is '
  '0.5–2.7 % lower up to the experimental gradients and 26 % lower at 10⁵ K m⁻¹ (Al‑3Cu‑5Nb‑0.1Fe), the difference coming '
  'from the geometric function and not from the elasto-plastic terms.')
P('Consistently with Eq. '
  f'({EQS["gep"]}), the surface tension of the elasto-plastic interface is composed of an elastic and a plastic '
  'component. The elastic component, $\\sigma_0+s^e$, is limited by the yield condition; the permanent areal strain '
  'contributes the plastic component, $\\left|\\sigma_0\\right|E_A^p$; and the derivative of the dissipated heat with respect to '
  'the areal strain, $-\\partial q_{diss}/\\partial E_A=\\sigma_y^s C_A/(C_A+H)$, relieves the tension:')
EQ(r'\tau_{SL}^{ep}=\left(\sigma_0+s^e\right)+\left|\sigma_0\right|E_A^p+\frac{C_A}{C_A+H}\sigma_y^s=\tau_{SL}^{el}+\sigma_y^s\quad\left(H=0\right)', 'tauep')
P('After yielding, the elastic component remains constant at $\\sigma_0-\\sigma_y^s$, the plastic component carries the '
  'growth of the tension with the gradient, and the total tension follows the elastic one, lowered in magnitude only by '
  '$\\sigma_y^s$.')
P('For comparison, three alternative closures were also examined: the release of the permanent area through the '
  f'reference radius of Eq. ({EQS["ref"]}), $\\gamma_0(r_{{eq,p}}/r)^2+\\left|\\sigma_0\\right|\\left|E_A\\right|$, which lowers the '
  'second-order tensor by 6 % at the experimental gradient of Al‑3Cu‑5Nb‑0.1Fe and by up to 29 % at 8×10⁴ K m⁻¹; the '
  'purely elastic part, $\\gamma_0(r_{eq,p}/r)^2-s$, in which the energy of the permanent area is lost; and the '
  '$\\gamma_0(r_{eq,p}/r)^2-s+w_p$ closure, in which only the plastic work is added. In the last two, the surface energy '
  'almost stops growing after the yield point and the nucleation barrier decreases with the thermal gradient, giving a '
  f'monotonically decreasing nucleation rate after yielding. The stored closure of Eq. ({EQS["gep"]}) was therefore '
  'adopted. At the liquidus temperature the interface does not harden, and $H=0$ is used.')
P('The radial derivative of the elasto-plastic surface energy, needed in Eqs. '
  f'({EQS["GT"]})–({EQS["G1"]}), follows from Eqs. ({EQS["tangent"]}) and ({EQS["gep"]}):')
EQ(r'\frac{\partial\gamma_{SL}^{ep}}{\partial r}=-\frac{2\gamma_0 r_{eq}^2}{r^3}-\frac{C_A H}{C_A+H}\frac{\partial E_A}{\partial r}-C_A\frac{\partial E_A^p}{\partial r}-\sigma_y^s\frac{\partial\alpha}{\partial r},\qquad \frac{\partial E_A^p}{\partial r}=-\frac{\partial\alpha}{\partial r}=\frac{C_A}{C_A+H}\frac{\partial E_A}{\partial r}', 'dgep')
P('during yielding ($s<0$), with $\\partial E_A^p/\\partial r=\\partial\\alpha/\\partial r=0$ and $ds/dE_A=C_A$ in the elastic '
  'range. The first term is the geometric term, the second the consistent tangent of the excess stress, the third the '
  'stored energy of the permanent deformation and the fourth the dissipated heat. For $H=0$, '
  '$\\partial\\gamma_{SL}^{ep}/\\partial r=-2\\gamma_0 r_{eq}^2/r^3-\\left(\\left|\\sigma_0\\right|-\\sigma_y^s\\right)\\partial E_A/\\partial r$, i.e., the '
  f'elastic derivative of Eq. ({EQS["dgel"]}) with $\\left|\\sigma_0\\right|$ reduced by $\\sigma_y^s$.')

H2('3.6. Elastic restoration driven by the rate of density change')
P('A specimen deformed dynamically shows an elastic deformation, an elastic restoration and a permanent deformation. '
  'In the nucleus, the recoverable (elastic) part of the surface energy is restored by the phase transformation itself, '
  'which changes the density of the system and unloads the interface, whereas the permanent deformation and its stored '
  'energy remain. The rate of volumetric deformation of the transformation is the rate of density change:')
EQ(r'\dot\varepsilon_v=\frac{1}{\rho}\frac{d\rho}{dt}=\beta_s\dot f_s=\beta_s\frac{\partial f_s}{\partial T}\dot T\cong\beta_s\frac{\dot T}{\Delta T_{SL}},\qquad \beta_s=\frac{\rho_S-\rho_L}{\rho_S}', 'epsv')
P('where $\\dot T=\\nabla T_L V_L$ is the local cooling rate, $f_s$ the solid fraction and $\\Delta T_{SL}$ the solidification '
  'interval, evaluated as the mean integral of $\\dot T\\,t_{SL}$ over the measured range of positions $P$:')
EQ(r'\Delta T_{SL}=\frac{1}{P_2-P_1}\int_{P_1}^{P_2}\dot T(P)\,t_{SL}(P)\,dP', 'DTSL')
P('Let $f$ be the restored fraction of the elastic energy. The restoration is produced by the density change at the '
  'rate $\\dot\\varepsilon_v(1-f)$ and is lost, by the reloading of the interface by the capillary pressure, with a '
  'relaxation time $\\tau_r$. In the linear (Onsager) regime of non-equilibrium thermodynamics [15, 16]:')
EQ(r'\frac{df}{dt}=\dot\varepsilon_v\left(1-f\right)-\frac{f}{\tau_r}\quad\Rightarrow\quad f=\frac{\dot\varepsilon_v\tau_r}{1+\dot\varepsilon_v\tau_r}', 'f')
P('for the quasi-steady state, $df/dt\\approx 0$. Only the elastic energy, $\\left|s\\right|=C_A\\left|E_A^e\\right|\\le\\sigma_y^s$, is '
  'restorable; the surface energy of the restored interface and the corresponding second-order tensor are')
EQ(r'\gamma_{SL}^{R}=\gamma_{SL}^{ep}-f\left|s\right|,\qquad \Gamma_R^{2nd}=\Gamma_{ep}^{2nd}-f\frac{\left|s\right|}{\left|\Delta S_V\right|}\ge\Gamma_{ep}^{2nd}-\frac{\sigma_y^s}{\left|\Delta S_V\right|}', 'gR')
P('Since the geometric term, $\\gamma_0(r_{eq}/r)^2>\\gamma_0$, and the stored energy of the permanent deformation are '
  'not restored, the thermal-field tensor of the plastically deformed interface is always larger than its equilibrium '
  'value, $\\gamma_0/\\left|\\Delta S_V\\right|$, and grows with the thermal gradient in the same way as the loaded tensor. The '
  'excess over equilibrium is the energy stored by the permanent surface deformation, which must be available for '
  'deformation-induced phase transformation, i.e., for the nucleation of new nuclei of another phase on the deformed '
  'interface. The restoration is negligible near equilibrium, $\\dot\\varepsilon_v\\tau_r\\ll 1$, and complete for fast density '
  'changes, $\\dot\\varepsilon_v\\tau_r\\gg 1$, but it is bounded by $\\sigma_y^s/\\left|\\Delta S_V\\right|$. Equation '
  f'({EQS["G1"]}) shows that the first-order tensor is lower than $\\Gamma^{{2nd}}$ by the logarithmic derivative of the '
  'surface energy; it is an approximation of the full quadratic solution, and its lower value is not a restoration of the '
  'stored energy.')

H2('3.7. Summary of the corrections of the surface energy')
TABLE('Table 1. Corrections of the surface energy of the nucleus introduced by the elasto-plastic surface.',
      ['Correction', 'Expression', 'Effect on the thermal-field tensor'],
      [['Areal yield of the excess stress', f'Eqs. ({EQS["yield"]})–({EQS["rm2"]})', 'End of the elastic range; limits $s$ and $\\tau$'],
       ['Permanent areal strain', f'Eq. ({EQS["split"]})', 'Stored in the interface; $r_{eq}$ unchanged'],
       ['Stored energy of the permanent deformation', f'Eq. ({EQS["gep"]})', '$\\left|\\sigma_0\\right|\\left|E_A^p\\right|$ kept in $\\gamma_{SL}$'],
       ['Hardening energy and dissipation', f'Eqs. ({EQS["psi"]}), ({EQS["wp"]}), ({EQS["gep"]})', '$\\psi_p$ stored; $q_{diss}$ leaves as heat'],
       ['Consistent tangent', f'Eqs. ({EQS["tangent"]}), ({EQS["dgep"]})', '$\\partial\\gamma_{SL}/\\partial r$ in $\\Gamma^{1st}$ and $\\Gamma^{2nd}$'],
       ['Elastic restoration', f'Eqs. ({EQS["epsv"]})–({EQS["gR"]})', 'Only $\\left|s\\right|\\le\\sigma_y^s$; $\\Gamma>\\gamma_0/\\left|\\Delta S_V\\right|$']],
      [5.4, 4.0, 6.6])

H2('3.8. Nucleation rate')
P('The nucleation rate of the formulation is [1]')
EQ(r'I=\frac{D}{a^2}\,\frac{4\pi r_C^2}{a^2}\,N_V\exp\left(\frac{\Delta G_C}{\Delta G_{C,ref}}\right)', 'I')
P(f'where $D$ is the self-diffusion coefficient, $a$ the jump distance, $N_V$ the density of modes of the new phase, '
  f'$\\Delta G_C$ the critical free energy of Eq. ({EQS["dG"]}) and $\\Delta G_{{C,ref}}$ its value at the first point of the '
  'sweep. The elasto-plastic surface enters the rate through $r_C$, $\\Delta T$ and $\\Delta G_C$.')

# ------------------------------------------------------------------ 4
H1('4. Methods')
P('The elastic solution of Refs. [1–3] for the FCC_A1 phase of the Al‑3Cu‑5Nb‑0.1Fe, Al‑0.8Si‑0.6Mg‑0.2Fe and '
  'Al‑3Cu‑0.5Mg‑0.2Fe alloys was kept unchanged, and the elastic surface of the main loop was replaced by the '
  f'elasto-plastic surface of Section 3. The radius was swept from $r_{{eq}}$ downwards; at each radius Eqs. '
  f'({EQS["rm1"]})–({EQS["rm2"]}) give $s$, $E_A^p$ and $\\alpha$, Eqs. ({EQS["gep"]}) and ({EQS["dgep"]}) give $\\gamma_{{SL}}$ and '
  f'$\\partial\\gamma_{{SL}}/\\partial r$, and Eq. ({EQS["DT"]}) gives $\\Delta T$ and $\\nabla T$ in closed form, for the elastic and '
  'the elasto-plastic surface alike. The sweep was extended until $\\Delta T=30$ K, which covers the whole range of '
  'solidification gradients up to about 10⁶ K m⁻¹, where the root of the second-order contact angle disappears. The material parameters are listed in Table 2; the areal yield stresses '
  'are of the order of a few percent of $\\left|\\sigma_0\\right|$, and $H=0$. Since $q_{diss}=\\sigma_y^s\\alpha$ vanishes both for '
  '$\\sigma_y^s\\to 0$ and for $\\sigma_y^s\\to\\infty$, the elasto-plastic tensor is bounded by the elastic one, and its deviation is '
  'largest at intermediate yield stresses.')
TABLE('Table 2. Parameters of the elastic and elasto-plastic Gurtin–Murdoch surface (FCC_A1). Common values: '
      '$\\lambda=50.34$ GPa, $\\mu=25.94$ GPa, $\\lambda_0=-5$ N m⁻¹, $\\mu_0=1.6$ N m⁻¹, $r_{eq}=6.3$ μm, $H=0$.',
      ['Alloy', '$T_L$ [K]', '$\\gamma_0$ [J m⁻²]', '$\\left|\\sigma_0\\right|$ [N m⁻¹]', '$\\sigma_y^s$ [N m⁻¹]', '$\\Delta S_V$ [J m⁻³ K⁻¹]', '$\\overline{\\nabla T}_{exp}$ [K m⁻¹]'],
      [['Al‑3Cu‑5Nb‑0.1Fe', '936.58', '0.296', '1.3761', '0.0386', '−9.799×10⁵', '3294'],
       ['Al‑0.8Si‑0.6Mg‑0.2Fe', '925.32', '0.183', '1.090', '0.1572', '−9.261×10⁵', '7902'],
       ['Al‑3Cu‑0.5Mg‑0.2Fe', '922.61', '0.169', '1.209', '0.0423', '−8.347×10⁵', '3708']],
      [3.6, 1.6, 1.9, 2.0, 2.0, 2.6, 2.3])
P('The implementation was verified in three ways: (a) the analytical derivative of Eq. '
  f'({EQS["dgep"]}) agrees with the finite-difference derivative of Eq. ({EQS["gep"]}) to a relative deviation of 3×10⁻⁶ away '
  'from the yield point; (b) with $\\sigma_y^s\\to\\infty$, every column of the elastic solution is reproduced exactly; and (c) '
  f'without the dissipated heat, Eq. ({EQS["gep"]}) coincides with the elastic surface energy to 2.5×10⁻¹⁶. For the elastic '
  'restoration (Al‑3Cu‑5Nb‑0.1Fe, the alloy with the local gradient measured along the casting), $\\nabla T_L(P)$ is the '
  'quadratic fit of the thermal data between 5 and 90 mm, $V_L=3.4P^{-0.56}$ mm s⁻¹ and $t_{SL}=4.93116P^{0.87606}$ s [17].')

# ------------------------------------------------------------------ 5
H1('5. Results and Discussion')
H2('5.1. End of the elastic range')
P('Table 3 lists the yield point of the interface, where the elastic range ends and the plastic range begins. The excess '
  'stress reaches the areal yield stress when the radius has decreased by only 1.5–7 % from $r_{eq}$, at undercoolings of '
  '0.09–0.14 K and thermal gradients of 560–925 K m⁻¹, below the experimental gradients of the three alloys. The nuclei '
  'that control the solidification of these alloys are therefore in the elasto-plastic range.')
TABLE('Table 3. Yield point of the elasto-plastic Gurtin–Murdoch surface (end of the elastic range) and second-order tensor at 10⁵ K m⁻¹.',
      ['Alloy', '$r/r_{eq}$ at yield', '$\\Delta T$ at yield [K]', '$\\nabla T$ at yield [K m⁻¹]', '$\\Gamma^{2nd}$ at 10⁵ K m⁻¹, elastic / elasto-plastic [m K]'],
      [['Al‑3Cu‑5Nb‑0.1Fe', '0.985', '0.114', '731', '5.693×10⁻⁶ / 5.644×10⁻⁶'],
       ['Al‑0.8Si‑0.6Mg‑0.2Fe', '0.930', '0.136', '925', '4.805×10⁻⁶ / 4.581×10⁻⁶'],
       ['Al‑3Cu‑0.5Mg‑0.2Fe', '0.980', '0.087', '562', '5.210×10⁻⁶ / 5.142×10⁻⁶']],
      [3.6, 2.2, 2.4, 2.6, 5.2])
P('With the closed form of Eq. '
  f'({EQS["DT"]}), both surfaces were solved beyond 10⁵ K m⁻¹, up to about 1.1×10⁶ K m⁻¹, where the root of the '
  'second-order contact angle disappears for both. At 10⁵ K m⁻¹ (ΔT ≈ 5 K), the elasto-plastic tensor is 0.9–4.7 % lower '
  'than the elastic one (Table 3).')

H2('5.2. Surface tension, surface energy and thermal-field tensor')
P('Figure 1 compares, curve by curve, the elastic and elasto-plastic solutions of Al‑3Cu‑5Nb‑0.1Fe as functions of the '
  'thermal gradient; the shaded band and the dash-dotted line mark the elastic range and its end. Table 4 lists the '
  'values at the experimental gradients of the three alloys. In the elastic range both solutions coincide. After '
  'yielding, the elastic component of the surface tension remains constant at $\\sigma_0-\\sigma_y^s$ (−1.415, −1.247 and '
  '−1.251 N m⁻¹), while its plastic component grows with the gradient, so that the total surface tension, Eq. '
  f'({EQS["tauep"]}), follows the elastic one (1.6, 5.3 and 1.7 % lower in magnitude at the experimental gradients). '
  'The surface energy, its derivative, '
  'the second-order tensor, the undercooling, the critical radius and the critical free energy of the elasto-plastic '
  'interface practically coincide with the elastic ones and grow with the gradient in the same way, because the '
  'permanent deformation is stored and only the dissipated heat is removed.')
FIG(os.path.join(HERE, 'Fig1_curves.png'),
    'Figure 1. Elastic (solid) and elasto-plastic (dashed) Gurtin–Murdoch solutions of the FCC_A1 phase of '
    'Al‑3Cu‑5Nb‑0.1Fe against the thermal gradient: surface tension (elastic + plastic), surface energy and its radial derivative, '
    'second-order thermal-field tensor, undercooling, critical radius, critical free energy, and homogeneous and '
    'heterogeneous nucleation rates. Shaded band: elastic range; dash-dotted line: yield point (end of the elastic '
    'range, 731 K m⁻¹); dotted line: mean experimental gradient.')
TABLE('Table 4. Elastic and elasto-plastic solutions at the mean experimental thermal gradient (elastic / elasto-plastic).',
      ['Quantity', 'Al‑3Cu‑5Nb‑0.1Fe (3294 K m⁻¹)', 'Al‑0.8Si‑0.6Mg‑0.2Fe (7902 K m⁻¹)', 'Al‑3Cu‑0.5Mg‑0.2Fe (3708 K m⁻¹)'],
      [['$\\tau_{SL}$, elastic + plastic [N m⁻¹]', '−1.968 / −1.937', '−2.090 / −1.980', '−1.840 / −1.808'],
       ['$\\gamma_{SL}$ [J m⁻²]', '1.047 / 1.041', '1.458 / 1.397', '0.916 / 0.908'],
       ['$\\partial\\gamma_{SL}/\\partial r$ [J m⁻³]', '−7.21×10⁵ / −7.09×10⁵', '−7.78×10⁵ / −7.24×10⁵', '−6.16×10⁵ / −6.02×10⁵'],
       ['$\\Gamma^{1st}$ [m K]', '3.887×10⁻⁷ / 3.901×10⁻⁷', '7.637×10⁻⁷ / 7.503×10⁻⁷', '4.171×10⁻⁷ / 4.182×10⁻⁷'],
       ['$\\Gamma^{2nd}$ [m K]', '1.069×10⁻⁶ / 1.063×10⁻⁶', '1.575×10⁻⁶ / 1.508×10⁻⁶', '1.097×10⁻⁶ / 1.088×10⁻⁶'],
       ['$\\Delta T$ [K]', '0.421 / 0.419', '0.791 / 0.774', '0.452 / 0.450'],
       ['$r_C$ [μm]', '5.081 / 5.067', '3.982 / 3.897', '4.853 / 4.832'],
       ['$\\Delta G_C$ [J]', '1.132×10⁻¹⁰ / 1.120×10⁻¹⁰', '9.68×10⁻¹¹ / 8.89×10⁻¹¹', '9.04×10⁻¹¹ / 8.88×10⁻¹¹']],
      [3.4, 4.2, 4.2, 4.2])
P('At the experimental gradients, the elasto-plastic $\\Gamma^{2nd}$ is 0.5, 4.2 and 0.9 % lower than the elastic one. The '
  'difference is the dissipated heat, $q_{diss}/\\left|\\Delta S_V\\right|$, and is largest for Al‑0.8Si‑0.6Mg‑0.2Fe, which has the '
  'largest areal yield stress and the largest experimental gradient. The first-order tensor changes by less than 2 %. '
  'The ratio $\\Gamma^{1st}/\\Gamma^{2nd}$ is 0.37, 0.50 and 0.38, i.e., the first-order approximation removes about half of '
  f'the second-order tensor through the derivative of the surface energy, Eq. ({EQS["G1"]}).')

H2('5.3. Nucleation rate')
P('Figure 2 shows the homogeneous nucleation rate of Eq. '
  f'({EQS["I"]}) against the thermal gradient. With the permanent deformation stored in the surface energy, the rate '
  'coincides with the elastic one in the elastic range and keeps its shape after yielding: it rises to a maximum and '
  'decreases. The maximum is 4, 19 and 7 % lower and lies at gradients 2, 19 and 6 % lower (Table 5); the largest '
  'difference again corresponds to the largest areal yield stress. The elasto-plastic surface therefore displaces the '
  'nucleation rate slightly, without changing its physical behaviour.')
FIG(os.path.join(HERE, 'Fig2_nucleation_rate.png'),
    'Figure 2. Homogeneous nucleation rate of the FCC_A1 phase against the thermal gradient, with the elastic and the '
    'elasto-plastic Gurtin–Murdoch surface. Shaded band: elastic range; dash-dotted lines: yield point; dotted lines: '
    'mean experimental gradients.')
TABLE('Table 5. Maximum of the homogeneous nucleation rate.',
      ['Alloy', 'Elastic: $I_{max}$ [m⁻³ s⁻¹] at $\\nabla T$ ($\\Delta T$)', 'Elasto-plastic: $I_{max}$ [m⁻³ s⁻¹] at $\\nabla T$ ($\\Delta T$)'],
      [['Al‑3Cu‑5Nb‑0.1Fe', '6.62×10⁴⁹ at 4953 K m⁻¹ (0.58 K)', '6.37×10⁴⁹ at 4873 K m⁻¹ (0.57 K)'],
       ['Al‑0.8Si‑0.6Mg‑0.2Fe', '1.01×10⁵⁰ at 4465 K m⁻¹ (0.51 K)', '8.19×10⁴⁹ at 3632 K m⁻¹ (0.43 K)'],
       ['Al‑3Cu‑0.5Mg‑0.2Fe', '1.43×10⁵⁰ at 5728 K m⁻¹ (0.64 K)', '1.33×10⁵⁰ at 5413 K m⁻¹ (0.61 K)']],
      [3.8, 6.1, 6.1])

H2('5.4. The whole range of solidification gradients')
P('Solidification experiments reach high thermal gradients, and the solution was therefore computed over the whole '
  'range, from the elastic range up to about 10⁶ K m⁻¹ (Figures 3–5, Table 6). Over this range, the second-order '
  'tensor and the surface energy of the elasto-plastic interface follow the elastic ones; the relative difference '
  'remains below 1 % for Al‑3Cu‑5Nb‑0.1Fe and Al‑3Cu‑0.5Mg‑0.2Fe and below 5 % for Al‑0.8Si‑0.6Mg‑0.2Fe. At 10⁵ K m⁻¹, '
  'where $\\Delta T\\approx 5$ K and $r_C\\approx 2$ μm, the surface energy reaches 4.2–5.5 J m⁻², i.e., 19–25 times its '
  'equilibrium value, and so does the thermal-field tensor relative to $\\gamma_0/\\left|\\Delta S_V\\right|$.')
P('The surface tension (Figure 5) shows how the two components share the load: after yielding, the elastic component '
  'is limited to $\\sigma_0-\\sigma_y^s$, and the plastic component, $\\left|\\sigma_0\\right|E_A^p$, carries the growth of the tension '
  'with the gradient (−2.4 to −3.0 N m⁻¹ at 10⁵ K m⁻¹). The total remains within 0.6–2.9 % of the elastic tension at '
  '10⁵ K m⁻¹. For Al‑0.8Si‑0.6Mg‑0.2Fe, which has the largest areal yield stress, the total tension shows a step of '
  f'$\\sigma_y^s$ at the yield point, the onset of the dissipation term of Eq. ({EQS["tauep"]}).')
FIG(os.path.join(HERE, 'Fig4_gibbs_thomson_range.png'),
    'Figure 3. Second-order thermal-field tensor of the FCC_A1 phase against the thermal gradient over the whole range, '
    'with the elastic and the elasto-plastic Gurtin–Murdoch surface, for the three alloys. Shaded band: elastic range; '
    'dash-dotted line: yield point; dotted lines: mean experimental gradient and equilibrium value '
    '$\\gamma_0/\\left|\\Delta S_V\\right|$; red: values at 10⁵ K m⁻¹.')
FIG(os.path.join(HERE, 'Fig5_surface_energy_range.png'),
    'Figure 4. Surface energy of the FCC_A1 nucleus against the thermal gradient over the whole range, elastic and '
    'elasto-plastic Gurtin–Murdoch surface, for the three alloys. Shaded band: elastic range; dash-dotted line: yield '
    'point; dotted lines: mean experimental gradient and equilibrium value $\\gamma_0$; red: values at 10⁵ K m⁻¹.')
FIG(os.path.join(HERE, 'Fig6_surface_tension_range.png'),
    'Figure 5. Surface tension of the FCC_A1 nucleus against the thermal gradient over the whole range, for the three '
    'alloys: elastic Gurtin–Murdoch surface and elasto-plastic surface, total (elastic + plastic) and its elastic, '
    '$\\sigma_0+s^e$, and plastic, $\\left|\\sigma_0\\right|E_A^p$, components. Shaded band: elastic range; dash-dotted line: yield '
    'point; red: values at 10⁵ K m⁻¹.')
TABLE('Table 6. Elastic and elasto-plastic solutions at 10⁵ K m⁻¹ (elastic / elasto-plastic).',
      ['Alloy', '$\\Gamma^{2nd}$ [m K]', '$\\gamma_{SL}$ [J m⁻²]', '$\\tau_{SL}$ [N m⁻¹]', '$\\Delta T$ [K]', '$r_C$ [μm]'],
      [['Al‑3Cu‑5Nb‑0.1Fe', '5.693×10⁻⁶ / 5.644×10⁻⁶', '5.579 / 5.531', '−4.362 / −4.336', '5.35 / 5.33', '2.13 / 2.12'],
       ['Al‑0.8Si‑0.6Mg‑0.2Fe', '4.805×10⁻⁶ / 4.581×10⁻⁶', '4.450 / 4.242', '−3.640 / −3.535', '4.91 / 4.80', '1.96 / 1.91'],
       ['Al‑3Cu‑0.5Mg‑0.2Fe', '5.210×10⁻⁶ / 5.142×10⁻⁶', '4.349 / 4.292', '−3.940 / −3.914', '5.12 / 5.08', '2.04 / 2.02']],
      [3.4, 3.6, 2.4, 2.6, 1.8, 1.8])

H2('5.5. Elastic restoration along the casting')
P('Figure 6 and Table 7 show the rate of density change and the elastic restoration of Section 3.6 along the horizontal '
  f'casting of Al‑3Cu‑5Nb‑0.1Fe, with $\\beta_s=0.0767$ ($\\rho_S=2737.04$, $\\rho_L=2527.05$ kg m⁻³) and $\\Delta T_{{SL}}=185.3$ K '
  f'from Eq. ({EQS["DTSL"]}). The rate of density change decreases by a factor of 55 from the cooled surface (5 mm) to '
  '90 mm. Since the whole interface has yielded, the restorable elastic energy is $\\left|s\\right|=\\sigma_y^s$ and the largest '
  'restoration of the tensor is $\\sigma_y^s/\\left|\\Delta S_V\\right|=3.94\\times 10^{-8}$ m K, 2–10 % of $\\Gamma_{ep}^{2nd}$. '
  'Even with complete restoration ($f=1$), the tensor of the deformed interface remains 1.2–6.3 times the equilibrium '
  'value, $\\gamma_0/\\left|\\Delta S_V\\right|=3.02\\times 10^{-7}$ m K, and decreases along the casting following the local gradient. '
  'The relaxation time $\\tau_r$, which sets $f$ between these bounds, and the analysis of the secondary dendrite arm '
  'spacing will be reported separately.')
FIG(os.path.join(HERE, 'Fig3_restoration.png'),
    'Figure 6. Elastic restoration along the casting of Al‑3Cu‑5Nb‑0.1Fe: (left) rate of density change; (right) '
    'second-order thermal-field tensor of the loaded elasto-plastic interface ($f=0$), range of the elastic restoration '
    '($0\\leq f\\leq 1$) and equilibrium value.')
TABLE('Table 7. Rate of density change and bounds of the elastic restoration along the casting of Al‑3Cu‑5Nb‑0.1Fe.',
      ['P [mm]', '$\\nabla T_L$ [K m⁻¹]', '$\\dot T$ [K s⁻¹]', '$\\dot\\varepsilon_v$ [s⁻¹]', '$\\Gamma_{ep}^{2nd}$, $f=0$ [m K]', '$\\Gamma_R^{2nd}$, $f=1$ [m K]', '$\\Gamma_R^{2nd}\\left|\\Delta S_V\\right|/\\gamma_0$'],
      [['5', '9431', '13.02', '5.39×10⁻³', '1.941×10⁻⁶', '1.902×10⁻⁶', '6.30'],
       ['20', '6820', '4.33', '1.79×10⁻³', '1.633×10⁻⁶', '1.593×10⁻⁶', '5.27'],
       ['40', '4072', '1.75', '7.27×10⁻⁴', '1.212×10⁻⁶', '1.173×10⁻⁶', '3.88'],
       ['60', '2161', '0.74', '3.07×10⁻⁴', '8.05×10⁻⁷', '7.65×10⁻⁷', '2.53'],
       ['80', '1088', '0.32', '1.32×10⁻⁴', '4.87×10⁻⁷', '4.47×10⁻⁷', '1.48'],
       ['90', '866', '0.24', '9.81×10⁻⁵', '4.06×10⁻⁷', '3.67×10⁻⁷', '1.21']],
      [1.3, 2.0, 1.7, 2.2, 2.8, 2.8, 2.5])

H2('5.6. Thermodynamic interpretation')
P('The elasto-plastic surface separates three energies of the interface that the elastic solution treats as one: the '
  'recoverable elastic energy, which is restored when the transformation unloads the interface; the energy of the '
  'permanent deformation, which is stored and remains available for nucleation; and the dissipated plastic work, which '
  'leaves the interface as heat and belongs to the path-dependent part of the thermal field. The Clausius–Duhem '
  f'inequality is satisfied identically by the flow rule of Eq. ({EQS["flow"]}) ($\\mathcal{{D}}^s=\\sigma_y^s\\dot\\lambda$), and the '
  f'restoration kinetics of Eq. ({EQS["f"]}) is a linear relaxation, so that both irreversible processes produce entropy. In '
  'the elastic range (low gradients), the elastic solution of Refs. [1–3] is recovered exactly. In the plastic range, the '
  'interface yields and stores its permanent deformation, so that the thermal-field tensor follows the elastic one, '
  'lowered only by the dissipated heat, and remains above its equilibrium value. This excess of the deformed interface '
  'over equilibrium is what makes possible the deformation-induced phase transformation, i.e., the nucleation of new '
  'nuclei of another phase on the plastically deformed interface, once the thermodynamic conditions are favourable.')

# ------------------------------------------------------------------------ 5.7 recrystallization
EP6 = os.path.join(HERE, '..', 'recristalizacao')
E64_AREA = 32.36e-10 ** 2      # cross-section of the 2048-atom cell [m^2] (the y-normal shear plane)


def e64_unload():
    """Each descarga_<pot>_<eps>.log in ../recristalizacao/ is one unloading run: it reverses the
    shear until the stress crosses zero. The state at the crossing has the elastic part removed:
    the residual strain and the stored energy per area there are the MD analogue of sigma_y*alpha.
    Returns the rows of Table 8."""
    import glob as _glob
    rows = []
    for log in sorted(_glob.glob(os.path.join(EP6, 'descarga_*.log'))):
        m0 = re.match(r'(?:\.\.[\\/]+)?recristalizacao[\\/]+descarga_(jelinek|mendelev)_([0-9.]+)\.log',
                      log.replace(os.sep, '/')) or \
             __import__('re').match(r'descarga_(jelinek|mendelev)_([0-9.]+)\.log', os.path.basename(log))
        if not m0:
            continue
        pot, eps0 = m0.group(1), float(m0.group(2))
        best = None
        for m in re.finditer(r'^\s*Step\s+Temp\s+Press\s+PotEng[^\n]*\n(.*?)(?:\n\s*\n|\n\s*Loop time|\Z)',
                             open(log, encoding='utf-8', errors='replace').read(), re.M | re.S):
            rr = [r.split() for r in m.group(1).splitlines()
                  if re.match(r'^\s*\d+\s', r) and len(r.split()) == 6]
            if rr and (best is None or len(rr) > len(best)):
                best = np.array([[float(x) for x in r] for r in rr])
        if best is None or np.ptp(best[:, 4]) == 0:
            continue
        step, pe, pxy = best[:, 0], best[:, 3] / 2048, best[:, 4]
        i0 = np.where((pxy[:-1] < 0) & (pxy[1:] >= 0))[0]
        if len(i0) == 0:
            continue
        i = int(i0[0])
        tc = float(np.interp(0.0, [pxy[i], pxy[i + 1]], [step[i], step[i + 1]]))
        pe_cross = float(np.interp(tc, [step[i], step[i + 1]], [pe[i], pe[i + 1]]))
        eps_cross = eps0 - 0.002 * tc * 0.002
        ref_m = re.search(r'^REF\s+(-[\d.]+)', open(os.path.join(EP6, f'deform_{pot}.log'),
                                                    encoding='utf-8', errors='replace').read(), re.M)
        ref = float(ref_m.group(1)) if ref_m else (-3.37336 if pot == 'mendelev' else -3.32191)
        stored = (pe_cross - ref) * 2048 * 1.602177e-19 / E64_AREA
        rows.append([pot, f'{eps0:.2f}', f'{eps_cross:.3f}', f'{pe_cross:.4f}', f'{stored:.3f}',
                     f'{0.157175 * (eps0 - 0.144):.4f}' if pot == 'jelinek' else '—'])
    return rows
    return rows


H2('5.7. Recrystallization driven by the stored plastic deformation: molecular-dynamics test')
P('The stored energy of the elasto-plastic interface is a thermodynamic state of the surface; it must therefore '
  'reappear as a driving potential in a second transformation. The natural test is the release: the plastically '
  'deformed solid is heated until the atomic diffusion activates and the recrystallization starts, consuming the '
  'stored deformation. The recrystallized grain nucleates in the deformed solid with the grain-boundary energy '
  '$\\gamma_{gb}$ in the place of the solid-liquid one, and the critical condition has the same quadratic structure '
  'of Section 2, now with the stored energy per unit volume, $\\sigma_y^s\\alpha/\\delta$, in the place of the bulk '
  'driving force $\\Delta S_V\\Delta T$:')
EQ(r'\Delta G_{rec}(r)=4\pi r^2\,\gamma_{gb}\,\frac{f(\theta)}{4}-\frac{4}{3}\pi r^3\,\frac{\sigma_y^s\,\alpha}{\delta},\qquad r^*=\frac{2\,\gamma_{gb}}{\sigma_y^s\,\alpha/\delta}\;\;\left(\theta=\pi\right)', 'grec')
P('Here $\\delta$ is the thickness of the deformed zone, which converts the areal storage of the interface into a '
  'volumetric driving force; it is measured in the atom maps of Figure 9 as the thickness of the high-energy band. '
  'Closed with this potential, the framework is self-consistent with the rules adopted in Section 3: the elastic '
  'part of the surface energy is restorable on its own (the elastic unloading), and the permanent part only leaves '
  'the interface consumed by the recrystallization, which returns the interface to equilibrium.')
P('The molecular-dynamics test (LAMMPS) runs FCC periodic cells of 2 048 atoms of Al‑0.8Si‑0.6Mg‑0.2Fe with the '
  'MEAM interatomic potential of Jelinek et al. [18] (the Cu of the parameterization is not included) and of pure '
  'Al with the EAM of Mendelev et al. [19] as the control. The cells are sheared at 300 K in strain steps of 10 ps '
  'at rates of 4–10 × 10⁹ s⁻¹, with a relaxation of 20 ps after each step; the stored energy is the potential '
  'energy per atom of the relaxed state relative to the undeformed reference, and $E_A^p=\\sigma_y^s/\\left|\\sigma_0\\right|=0.144$ '
  'of the model (the yield strain of the Al-Si-Mg-Fe interface) marks Figure 7.')
FIG(os.path.join(HERE, 'Fig7_E6_stored.png'),
    'Figure 7. The molecular-dynamics test of the stored plastic deformation: (left) stored energy per atom against '
    'the shear strain, with the yield strain of the model, $E_A^p=\\sigma_y^s/\\left|\\sigma_0\\right|=0.144$, marked; (right) the '
    'shear stress against the strain, with the flow. Data: Jelinek MEAM of the Al-0.8Si-0.6Mg-0.2Fe alloy and '
    'Mendelev EAM of pure Al as the control.')
P('The flow appears exactly where the elasto-plastic surface predicts the storage. The control yields at the model '
  'strain: the stored energy peaks at 0.046 eV/atom at $\\varepsilon=0.14$ and some 75 % of it is released by the flow of the '
  'crystal, whose shear stress falls from 5.9 to 1.9 GPa. The solute-strengthened alloy yields later, between '
  '0.20 and 0.30: its stored energy peaks at 0.051 eV/atom at 0.20, with a shear stress that falls from 4.0 GPa to '
  '0.6 GPa across the flow. What remains of the stored energy after the flow is the permanent deformation that the '
  'elasto-plastic interface claims to be available for nucleation.')
FIG(os.path.join(HERE, 'Fig8_E6_anneal.png'),
    'Figure 8. Annealing steps of the deformed cells: potential energy along each temperature step (ramp 300 K to T, '
    '100 ps at the plate). Left: the alloy, with the violent threshold between 700 and 800 K. Right: the pure-Al '
    'control, with the recovery at 800-900 K (the release of the stored energy, about 0.05 eV per atom).')
P('The release by heating (Figure 8) shows the stored potential leaving the deformed solid: at 800–900 K the '
  'control recovers about 0.05 eV per atom (about 5 × 10⁸ J m⁻³, three orders above the classic cold-work values - '
  'consistent with the violent nanoscale shear and with the threshold) inside 100 ps, while the alloy crosses a violent '
  'threshold between 700 and 800 K, where the deformed state decomposes (a catastrophic release of the stored '
  'energy, with atoms ejected at the very front in the small periodic cell). In 100 ps of molecular dynamics no '
  'perfect FCC grain nucleates yet - the recrystallization onset lives between the threshold of 800 K and the '
  'longer timescales - so the measurement of the release relies on the energy of the interface states rather than '
  'on grain shapes.')
FIG(os.path.join(HERE, 'Fig9_E6_field.png'),
    'Figure 9. The stored deformation as an energy field, atom by atom (energy per atom, colored: bright is stored '
    'in the defect network). From left to right: the sheared state of the control ($\\varepsilon=0.14$, the yield state of the '
    'model); the same cell after 100 ps at 500 K (still stored) and at 800 K (released); and the sheared state of the '
    'Al-Si-Mg-Fe alloy. The thickness of the bright band measures $\\delta$, the deformed zone of Eq. '
    f'({EQS["grec"]}).')
e64 = e64_unload()
if e64:
    TABLE('Table 8. The elastic unloading of the sheared cells: the energy remaining at zero shear stress '
          '(the stored permanent deformation, the MD analogue of $\\sigma_y^s\\alpha$) against the prediction of the '
          'stored_total mode for the alloy.',
          ['Surface potential', 'Base strain $\\varepsilon_0$', 'Strain after unloading',
           '$E$ at zero stress [eV/atom]', 'Stored remainder [J m⁻²]',
           'Model $\\sigma_y^s(\\varepsilon_0-E_A^p)$ [J m⁻²]'],
          [[r[0], r[1], r[2], r[3], r[4],
            '—' if r[0] == 'mendelev' else f'{0.157175 * (float(r[1]) - 0.144):.4f}'] for r in e64],
          [3.6, 2.4, 2.6, 2.6, 2.6, 2.8])
    P('The elastic unloading validates the measurement itself: the state at the peak of the stress (the elastic one, '
      '$\\varepsilon=0.14$) returns to zero stress with no stored energy left, below the equilibration noise (about 0.007 J m⁻² '
      'between the two machines). After the flow, the remainder is real and grows with the plastic history: 0.094 J m⁻² at '
      '$\\varepsilon=0.20$ and 0.184 J m⁻² at 0.30 in the pure-Al control, 1.3 J m⁻² in the alloy - 4 × 10⁸ J m⁻³ in the cell of '
      '32 Å, three orders of magnitude above the classic cold-work values of aluminium, the signature of the violent '
      'nanoscale shear. Per unit of remaining plastic strain the storage gives an effective areal flow stress of 1.7 N m⁻¹ '
      '(pure Al) and 10.3 N m⁻¹ (alloy): the solute-strengthened interface stores 6–13 times more energy for the same '
      'permanent deformation, which is what allows it to keep the larger thermal-field tensor of Section 5. The prediction '
      'of the stored_total mode with the $\\sigma_y^s$ fitted to the SDAS data (0.157 N m⁻¹) sits below these MD remainders of '
      'the deep states; the bridge between the two is the thickness $\\delta$ of the deformed zone of Eq. '
      f'({EQS["grec"]}), measured from the energy profiles: it converts the areal storage of the interface into the '
      'volumetric one the cell carries. This calibration of $\\sigma_y^s$ and $\\delta$ against the molecular-dynamics states '
      'is the natural continuation of the test.')
P('The size consistency of the measurement was tested against a 16 384-atom cell: the stored energy and the '
  'flow of the big cell reproduce the 2 048-atom ones within 3-5 % (stored 0.027, 0.049 and 0.040 eV/atom at '
  '0.14, 0.20 and 0.30 - the last one after the flow released the peak - shear stresses 3.49/3.52 GPa and '
  '4.10/3.97 GPa in the same strains), so the stored deformation of the sheared state is not a size effect at '
  'the scales tested.')


def g3_events():
    """The fine frontier: each recozer_jelinek_T<T>_s<S>.log is a 100 ps anneal of the deformed
    alloy cell with a different thermal seed. The catastrophic decomposition shows up as a jump
    of the potential energy and lost atoms. Returns the Table 9 rows (T, seed, event)."""
    import glob as _glob
    out = []
    for log in sorted(_glob.glob(os.path.join(EP6, 'gpu2', 'recozer_jelinek_T*_s*.log')) +
                      _glob.glob(os.path.join(EP6, 'recozer_jelinek_T*_s*.log'))):
        mm = re.match(r'recozer_jelinek_T([0-9.]+)_s([0-9]+)\.log', os.path.basename(log))
        if not mm:
            continue
        T, S = float(mm.group(1)), int(mm.group(2))
        if [(r[0], r[1]) for r in out] and (T, S) in [(float(r[0]), int(r[1])) for r in out]:
            continue
        best = None
        for m in re.finditer(r'^\s*Step\s+Temp\s+Press\s+PotEng[^\n]*\n(.*?)(?:\n\s*\n|\n\s*Loop time|\Z)',
                             open(log, encoding='utf-8', errors='replace').read(), re.M | re.S):
            rr = [r.split() for r in m.group(1).splitlines()
                  if re.match(r'^\s*\d+\s', r) and len(r.split()) == 6]
            if rr and (best is None or len(rr) > len(best)):
                best = np.array([[float(x) for x in r] for r in rr])
        if best is None:
            continue
        step, temp, pe = best[:, 0], best[:, 1], best[:, 3]
        dpe = np.diff(pe)
        i = int(np.argmax(dpe))
        raw = open(log, encoding='utf-8', errors='replace').read()
        lost = re.search(r'Lost atoms: original (\d+) current (\d+)', raw)
        n_lost = int(lost.group(1)) - int(lost.group(2)) if lost else 0
        if n_lost:
            out.append([f'{T:g}', str(S), f'{step[i] * 0.002:.0f}', f'{dpe[i]:+.0f}',
                        f'{temp.max():.0f}', str(n_lost)])
        else:
            out.append([f'{T:g}', str(S), 'limpo (100 ps)', '—', '—', '0'])
    return out


def long_run_row():
    """The 1 ns long anneal at 750 K (seed 777; its log is gpu2/longo_T750_5.log): the big cell the
    chain's shared filename substituted; the decay at 2 ps with 2 atoms lost."""
    f = os.path.join(EP6, 'gpu2', 'longo_T750_5.log')
    if not os.path.exists(f):
        return None
    raw = open(f, encoding='utf-8', errors='replace').read()
    lost = re.search(r'Lost atoms: original (\d+) current (\d+)', raw)
    if not lost:
        return None
    best = None
    for m in re.finditer(r'^\s*Step\s+Temp\s+Press\s+PotEng[^\n]*\n(.*?)(?:\n\s*\n|\n\s*Loop time|\Z)',
                         raw, re.M | re.S):
        rr = [r.split() for r in m.group(1).splitlines() if re.match(r'^\s*\d+\s', r) and len(r.split()) == 6]
        if rr and (best is None or len(rr) > len(best)):
            best = np.array([[float(x) for x in r] for r in rr])
    if best is None:
        return None
    step, temp, pe = best[:, 0], best[:, 1], best[:, 3]
    dpe = np.diff(pe)
    i = int(np.argmax(dpe))
    big = int(lost.group(1)) > 3000
    return ['750', '777 (longo, ' + ('célula grande' if big else 'célula pequena') + ')',
            f'{step[i] * 0.002:.0f}', f'{dpe[i]:+.0f}', f'{temp.max():.0f}',
            str(int(lost.group(1)) - int(lost.group(2)))]


H2('5.8. The stochastic onset: the fine frontier with thermal seeds')
P('The threshold of Figure 8 is not a state but a rate: at the plate of 750-800 K the deformed cell of the alloy '
  'decays catastrophically with a jump of the potential energy, a violent thermal spike and lost atoms, and the '
  'onset depends on the thermal history of the anneal. Table 9 gives the statistics of the fine frontier: 6 of the '
  '9 anneals decayed inside the 100 ps window, with onset times from 17 to 66 ps and no simple temperature trend '
  'at this 25 K spacing - the decay is probabilistic in the small cell, not a sharp threshold; the same '
  'temperature produced both clean and decayed cells among the seeds. The long step (1 ns, big cell) decayed '
  '2 ps into its ramp: in the 16 384-atom cell the decay at 750 K is almost immediate, which closes the map: '
  'the onset of the decay sits between 700 and 750 K.')
g3e = g3_events()
_long = long_run_row()
if _long:
    g3e.append(_long)
if g3e:
    TABLE('Table 9. The fine frontier of the decay of the deformed Al‑0.8Si‑0.6Mg‑0.2Fe cell '
          '(first seeds; the event = the jump of the potential energy, with lost atoms in the small cell).',
          ['T [K]', 'Seed', 'Event time [ps]', 'PE jump [eV]', 'T max [K]', 'Atoms lost'],
          g3e, [1.6, 1.3, 2.4, 2.4, 2.2, 2.2])
    FIG(os.path.join(HERE, 'Fig10_E6_onset.png'),
        'Figure 10. The onset times of the catastrophic decay across the fine frontier of the deformed cell '
        '(750-800 K, three thermal seeds): 6 of the 9 anneals decayed within 17-66 ps and the clean cells '
        'reached the full 100 ps window (squares) - no simple temperature trend at this spacing: the decay '
        'is probabilistic in the small cell. The big cell (16 384 atoms) decayed 2 ps into the anneal.')

# ------------------------------------------------------------------ 6
H1('6. Conclusions')
for c in [
    'The Gurtin–Murdoch surface was generalized to elasto-plasticity within non-equilibrium thermodynamics with internal '
    'variables: additive split of the areal strain, surface free energy with hardening, Clausius–Duhem inequality, areal '
    'yield condition, associative flow from maximum dissipation and return mapping with consistent tangent.',
    'The surface energy of the nucleus contains elastic and plastic terms: the permanent surface deformation complements '
    'the elastic one, its energy is stored in the interface and available for nucleation, and only the heat dissipated '
    'by the plastic flow leaves it. The elastic solution is recovered exactly in the elastic range.',
    'The interface of the FCC_A1 nuclei of the three alloys yields at 560–925 K m⁻¹, which marks the end of the elastic '
    'range and lies below the experimental gradients. Both surfaces are solved up to about 1.1×10⁶ K m⁻¹.',
    'The elasto-plastic thermal-field tensor grows with the thermal gradient as the elastic one: it is 0.5–4.2 % lower at '
    'the experimental gradients and 0.9–4.7 % lower at 10⁵ K m⁻¹, the difference being the dissipated heat. The '
    'nucleation rate keeps its shape, with a maximum 4–19 % lower.',
    'The surface tension of the deformed interface is composed of an elastic component, limited by the areal yield '
    'stress, and a plastic component of the permanent deformation; the total follows the elastic tension within '
    '0.6–5.3 %. Over the whole range of solidification gradients, up to 10⁶ K m⁻¹, the thermal-field tensor, the surface '
    'energy and the surface tension of the elasto-plastic interface follow the elastic ones; at 10⁵ K m⁻¹ the surface '
    'energy is 19–25 times its equilibrium value.',
    'Only the elastic part of the surface energy is restored by the rate of density change of the transformation, by at '
    'most $\\sigma_y^s/\\left|\\Delta S_V\\right|$. The thermal-field tensor of the plastically deformed interface therefore remains above '
    'its equilibrium value, as required for deformation-induced nucleation of new phases. The first-order tensor is an '
    'approximation of the full quadratic solution and does not represent a restoration of the stored energy.',
    'A molecular-dynamics test of the Al‑0.8Si‑0.6Mg‑0.2Fe alloy and of pure Al confirms the storage: the stored '
    'energy peaks at the yield strain of the model, the flow releases most of it, and the heating releases the '
    'remaining stored deformation as the driving potential of the recrystallization, which the same quadratic '
    'structure of Section 2 describes with $\\sigma_y^s\\alpha/\\delta$ in the place of $\\Delta S_V\\Delta T$. The deep '
    'state of the control (0.30) was unloaded independently on the two machines: 0.190 J m⁻² on the GPU box and '
    '0.184 J m⁻² on the laptop, an agreement of 4 %.',
]:
    par = doc.add_paragraph(style='List Number')
    par.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    rich(par, c)

H1('Declarations')
P('Competing interests. The author declares that there are no competing financial interests or personal relationships '
  'that may have influenced the work reported in this study.')
P('Authors’ contributions. Ferreira developed the formalism, derived the proposed equations, performed the computations, '
  'and wrote and revised this paper.')
P('Funding. The author acknowledges the financial support provided by CAPES (Coordenação de Aperfeiçoamento de Pessoal '
  'de Nível Superior – Brazil), Finance Code 001 and Grant 88881.707312/2022-01, CNPq (National Council for Scientific '
  'and Technological Development – Brazil), Grant 301502/2022-6, and FAPESPA (Fundação Amazônia de Amparo a Estudos e '
  'Pesquisas).')
P('Data availability. The Python implementation of the elasto-plastic surface and the data supporting the findings of '
  'this study are available in the GitHub repository https://github.com/ileaof/gm-elastoplastic-nucleation.')

H1('References')
refs = [
    'I.L. Ferreira. Assessment of thermodynamic variables affecting phase nucleation. Physica B: Condens. Matter 695 (2024) 416494.',
    'I.L. Ferreira. Thermal gradient affecting phase Gibbs free-energy. Materials (2026), in press.',
    'K.F. Kelton, A.L. Greer. Nucleation in Condensed Matter: Applications in Materials and Biology. Pergamon '
    'Materials Series, Vol. 15, Elsevier, Amsterdam, 2010.',
    'L.P. Marques, S.Q. Mesquita, A.C. Maciel, R.B. Costa, R.S. França, O.F.L. da Rocha, I.L. Ferreira. Effect of cooling rate on the secondary dendritic spacing of a horizontally solidified 6xxx series aluminum alloy. Mater. Res. 28 (2025) e20250045.',
    'I.L. Ferreira, A.L.S. Moreira. On the continuous mechanics first- and second-order formulations for nonequilibrium nucleation: derivation and applications. Int. J. Thermophys. 44 (2023) 72.',
    'M.E. Gurtin, A.I. Murdoch. A continuum theory of elastic material surfaces. Arch. Ration. Mech. Anal. 57 (1975) 291–323.',
    'M.E. Gurtin, A.I. Murdoch. Surface stress in solids. Int. J. Solids Struct. 14 (1978) 431–440.',
    'G.I. Taylor, H. Quinney. The latent energy remaining in a metal after cold working. Proc. R. Soc. Lond. A 143 (1934) 307–326.',
    'B.D. Coleman, W. Noll. The thermodynamics of elastic materials with heat conduction and viscosity. Arch. Ration. Mech. Anal. 13 (1963) 167–178.',
    'B.D. Coleman, M.E. Gurtin. Thermodynamics with internal state variables. J. Chem. Phys. 47 (1967) 597–613.',
    'J.R. Rice. Inelastic constitutive relations for solids: an internal-variable theory and its application to metal plasticity. J. Mech. Phys. Solids 19 (1971) 433–455.',
    'J.C. Simo, T.J.R. Hughes. Computational Inelasticity. Springer, New York, 1998.',
    'M.E. Gurtin, E. Fried, L. Anand. The Mechanics and Thermodynamics of Continua. Cambridge University Press, Cambridge, 2010.',
    'J. Lubliner. Plasticity Theory. Macmillan, New York, 1990.',
    'L. Onsager. Reciprocal relations in irreversible processes. I. Phys. Rev. 37 (1931) 405–426.',
    'S.R. de Groot, P. Mazur. Non-Equilibrium Thermodynamics. North-Holland, Amsterdam, 1962.',
    'Mendes et al. (2023) [complete reference to be inserted by the author: Al‑3Cu‑5Nb‑0.1Fe horizontal transient solidification].',
    'B. Jelinek et al. [complete reference to be confirmed by the author: modified embedded-atom method potential '
    'of the Al-Si-Mg-Cu-Fe alloy system]. Phys. Rev. B 85 (2012) 245102.',
    'A. M. Mendelev, M. J. Kramer, C. A. J. Becker, M. Asta. Analysis of semi-empirical interatomic potentials appropriate for '
    'simulation of crystalline and liquid Al and Cu. Philos. Mag. 88 (2008) 1723-1750.',
]
for i, rf in enumerate(refs, 1):
    par = doc.add_paragraph()
    par.paragraph_format.line_spacing = 1.15
    par.paragraph_format.left_indent = Cm(0.8)
    par.paragraph_format.first_line_indent = Cm(-0.8)
    r = par.add_run(f'[{i}] {rf}')
    r.font.size = Pt(11)

_SUP_MAP = {'⁰': '0', '⁴': '4', '⁵': '5', '⁶': '6', '⁷': '7',
            '⁹': '9', '⁻': '-'}
_NB_HYPHEN = '‑'          # non-breaking hyphen: absent from every journal XML entity set
_SUP_CHARS = ''.join(_SUP_MAP)  # ¹ ² ³ (Latin-1, entity-backed) are left as they are


def _sup_segments(text):
    """Split text into (is_superscript, characters) segments: the Unicode superscript
    digits/minus become the ASCII character; everything else is a normal segment."""
    segs, kind, buf = [], None, ''
    for ch in text:
        if ch in _SUP_MAP:
            k, c = True, _SUP_MAP[ch]
        else:
            k, c = False, ('-' if ch == _NB_HYPHEN else ch)
        if kind is None or k != kind:
            if buf:
                segs.append((kind, buf))
            kind, buf = k, c
        else:
            buf += c
    if buf:
        segs.append((kind, buf))
    return segs


def sanitize_text(doc):
    """ScholarOne rejects plain text carrying characters without a classic XML entity: the
    non-breaking hyphen U+2011 (it fills the alloy designations) and the Unicode superscript
    digits/minus (the 10^-6 K/m of the text and captions). This pass rewrites every plain-text
    run of the body (tables included): U+2011 -> '-', Unicode superscripts -> ASCII inside a
    real superscript run (w:vertAlign). Runs that mix tabs or breaks with text are left as they
    are (the builder never puts these characters there). Called once, right before the save."""
    from docx.text.run import Run
    space = '{http://www.w3.org/XML/1998/namespace}space'
    keep = {qn('w:lastRenderedPageBreak')}
    fixed = 0
    for t_el in list(doc.element.body.iter(qn('w:t'))):
        text = t_el.text or ''
        if _NB_HYPHEN not in text and not any(c in _SUP_CHARS for c in text):
            continue
        r_el = t_el.getparent()
        content = [c for c in r_el if c.tag != qn('w:rPr')]
        if r_el is None or r_el.tag != qn('w:r') or \
                any(c.tag not in keep | {qn('w:t')} for c in content):
            print('sanitize_text: skipped a run:', repr(text[:60]))
            continue
        segs = _sup_segments(text)
        if len(segs) == 1 and segs[0][0] is False:
            t_el.text, t_el.attrib[space] = segs[0][1], 'preserve'
            fixed += 1
            continue
        for j, (is_sup, seg) in enumerate(segs):
            nr = copy.deepcopy(r_el)
            if j > 0:                       # page-break markers stay in the first split run
                for marker in list(nr):
                    if marker.tag in keep:
                        nr.remove(marker)
            for other in nr.findall(qn('w:t'))[1:]:
                nr.remove(other)
            tt = nr.find(qn('w:t'))
            tt.text, tt.attrib[space] = seg, 'preserve'
            if is_sup:
                Run(nr, None).font.superscript = True
            else:
                rpr = nr.find(qn('w:rPr'))
                if rpr is not None:
                    for va in rpr.findall(qn('w:vertAlign')):
                        rpr.remove(va)
            r_el.addprevious(nr)
        r_el.getparent().remove(r_el)
        fixed += 1
    return fixed


n_sanitize = sanitize_text(doc)
doc.save(OUT)
print('saved', OUT, 'equations', eq_no[0], ', sanitized runs', n_sanitize)
