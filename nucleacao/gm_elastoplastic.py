"""
Elasto-plastic Gurtin-Murdoch (GM) surface for the nucleation models of I.L. Ferreira.

Same formulation as gm-elastoplastic-nucleation (docs/FORMULATION.md, sections 5-6):

    elastic GM (original scripts):  s = -2 sigma_0 ln(y0/y) = C_A E_A,  E_A = 2 ln(y/y0),  C_A = sigma_0
                                    gamma = gamma_0 (r_eq/r)^2 - s
    areal yield on the excess stress |s| <= sigma_y + H alpha (return mapping, r decreasing from r_eq)
    plastic area moves the equilibrium reference: r_eq,p = r_eq exp(E_A^p / 2)
    plastic work w_p = psi_p + q_diss,  psi_p = H alpha^2 / 2,  q_diss = sigma_y alpha

Settings (environment variables, defaults in brackets):
    GM_SURFACE_MODEL      elastic | elastoplastic                         [elastoplastic]
    GM_SIGMA_YA0          areal yield stress of the interface [N/m]       [value given by the script]
    GM_H_A                linear hardening [N/m]                          [0]
    GM_PLASTIC_REFERENCE  1: plastic area updates the reference, 0: fixed [0 for stored_total, 1 otherwise]
    GM_DISSIPATION_MODE   stored | thermo | barrier | restored | restored_barrier   [mode_default of the script]
                          barrier: gamma = gamma_0 (r_eq,p/r)^2 - s + w_p (elastic + plastic deformation energy)
                          stored_total: gamma = gamma_0 (r_eq/r)^2 + |sigma_0| (|E_A^e| + |E_A^p|) - q_diss, fixed
                          reference: the permanent deformation is stored, only the dissipated heat leaves (author)
                          complement: gamma = gamma_0 (r_eq,p/r)^2 - sigma_0 (E_A^e + E_A^p); the permanent
                          deformation complements the elastic one and its energy is stored for nucleation
    GM_DT_MAX             end of the r-sweep, undercooling [K]            [30, ~1e6 K/m]
GM_SURFACE_MODEL=elastic reproduces the original elastic scripts exactly.
"""
import os
from math import log, exp, pi


class GMElastoPlasticSurface:
    """
    surface_form: 'legacy'         gamma = gamma_0 (r_eq/r)^2 - s                     (gamma_func_r of most scripts)
                  'metric_2sphere' gamma = 3 gamma_0 / (1 + 2 r^2/r_eq^2) - s          (metric form of the *_dev_gam script)
    In both forms the plastic area replaces r_eq by r_eq,p in the gamma_0 term.
    """
    def __init__(self, req, gamma_0, sigma_0, lambda_, mu, lambda0, mu0, sigma_y_default=0.0695,
                 surface_form='legacy', mode_default='stored'):
        self.surface_form = surface_form
        self.mode_default = mode_default
        self.req, self.gamma_0, self.sigma_0 = req, gamma_0, sigma_0
        self.lambda_, self.mu, self.lambda0, self.mu0 = lambda_, mu, lambda0, mu0
        self.model = os.environ.get('GM_SURFACE_MODEL', 'elastoplastic')
        self.sigma_y = float(os.environ.get('GM_SIGMA_YA0', sigma_y_default))
        self.H = float(os.environ.get('GM_H_A', 0.0))
        self.mode = os.environ.get('GM_DISSIPATION_MODE', mode_default)
        # stored_total keeps the permanent area stored in the interface: the reference is not updated
        self.ref_update = os.environ.get('GM_PLASTIC_REFERENCE', '0' if self.mode == 'stored_total' else '1') == '1'
        self.DT_max = float(os.environ.get('GM_DT_MAX', 30.0))   # whole range of the sweep (~1e6 K/m)
        self.C_A = sigma_0
        self.EAp = self.alpha = self.wp = self.psi_p = self.q_diss = 0.0
        self.yielded_at = None
        # per step: (tau elastic component, plastic component |sigma_0| E_A^p, dissipation term, total tau)
        self.hist = []

    @property
    def tag(self):
        if self.model == 'elastic':
            return 'elastic'
        ref = 'refupd' if self.ref_update else 'refix'
        return f'EP_{self.mode}_{ref}_sy{self.sigma_y:g}_H{self.H:g}'

    def describe(self):
        if self.model == 'elastic':
            return 'elastic Gurtin-Murdoch surface (original)'
        return (f'elasto-plastic Gurtin-Murdoch surface: sigma_y = {self.sigma_y:g} N/m, H = {self.H:g} N/m, '
                f'reference {"updated by the plastic area" if self.ref_update else "fixed"}, energy mode {self.mode}')

    def _areal_strain(self, r):
        y = (3 * self.lambda_ + 2 * self.mu) * r + 2 * (self.lambda0 + self.mu0)
        y0 = (3 * self.lambda_ + 2 * self.mu) * self.req + 2 * (self.lambda0 + self.mu0)
        return 2 * log(y / y0)

    def _dEA_dr(self, r):
        # tangent of dgammadr_func_r of the original scripts
        alpha = 2 * (self.sigma_0 + 2 * (self.lambda0 + self.mu0)) / (r * (3 * self.lambda_ + 2 * self.mu))
        return 2 / (r * (1 + alpha))

    def step(self, r):
        """Advance the surface to radius r (r decreasing). Returns (tau, s, gamma, dgamma/dr)."""
        C, H = self.C_A, self.H
        E_A, dEA = self._areal_strain(r), self._dEA_dr(r)
        s = C * (E_A - self.EAp)
        k = C
        dalpha_dEA = 0.0
        yielding = False
        if self.model == 'elastoplastic':
            sy_n = self.sigma_y + H * self.alpha
            f_tr = abs(s) - sy_n
            if f_tr > 0.0:
                dlam = f_tr / (C + H)
                sgn = 1.0 if s > 0.0 else -1.0
                s -= C * dlam * sgn
                self.EAp += dlam * sgn
                self.alpha += dlam
                self.wp += (sy_n + 0.5 * H * dlam) * dlam
                self.q_diss += self.sigma_y * dlam
                self.psi_p = 0.5 * H * self.alpha ** 2
                k = C * H / (C + H)
                dalpha_dEA = sgn * C / (C + H)
                yielding = True
                if self.yielded_at is None:
                    self.yielded_at = r
        tau = -self.sigma_0 + s
        req_p = self.req * exp(self.EAp / 2) if (self.model == 'elastoplastic' and self.ref_update) else self.req
        dEAp = (C / (C + H)) * dEA if (yielding and self.ref_update) else 0.0
        if self.surface_form == 'metric_2sphere':
            x = (r / req_p) ** 2
            dx = 2 * r / req_p ** 2 - x * dEAp          # d(r^2/r_eq,p^2)/dr, with d r_eq,p/dr = r_eq,p dE_A^p/dr / 2
            gam = 3 * self.gamma_0 / (1 + 2 * x) - s
            dgam = -6 * self.gamma_0 / (1 + 2 * x) ** 2 * dx - k * dEA
        else:
            gam = self.gamma_0 * (req_p / r) ** 2 - s
            dgam = -2 * self.gamma_0 * req_p ** 2 / r ** 3 + self.gamma_0 * req_p ** 2 / r ** 2 * dEAp - k * dEA
        if self.model == 'elastoplastic':
            if self.mode == 'thermo':
                gam += self.psi_p
                dgam += H * self.alpha * dalpha_dEA * dEA
            elif self.mode == 'barrier':
                gam += self.wp
                dgam += (self.sigma_y + H * self.alpha) * dalpha_dEA * dEA
            elif self.mode == 'complement':
                # the permanent surface deformation complements the elastic one and its energy is stored:
                # the plastic area is also worked by the residual tension sigma_0 (as the elastic area in -s)
                #   gamma = gamma_0 (r_eq,p/r)^2 + sigma_0 (|E_A^e| + |E_A^p|) = gamma_0 (r_eq,p/r)^2 - sigma_0 E_A
                dEAp_r = (C / (C + H)) * dEA if yielding else 0.0
                gam += -C * self.EAp
                dgam += -C * dEAp_r
            elif self.mode == 'stored_total':
                # the permanent deformation is stored in the interface (energy of the residual tension on the
                # permanent area) and only the dissipated heat q_diss = sigma_y alpha leaves it:
                #   gamma = gamma_0 (r_eq/r)^2 + |sigma_0| (|E_A^e| + |E_A^p|) - q_diss
                dEAp_r = (C / (C + H)) * dEA if yielding else 0.0
                gam += -C * self.EAp - self.q_diss
                dgam += -C * dEAp_r - self.sigma_y * dalpha_dEA * dEA
            elif self.mode in ('restored', 'restored_barrier'):
                gam, dgam = self.gamma_0, 0.0
                if self.mode == 'restored_barrier':
                    gam += self.wp
                    dgam += (self.sigma_y + H * self.alpha) * dalpha_dEA * dEA
        # surface tension composed of the elastic and plastic components (author, 2026-10-23), consistent with
        # gamma = gamma_el - sigma_y alpha: tau_total = (sigma_0 + s^e) + |sigma_0| E_A^p + d(-q_diss)/dE_A
        if self.model == 'elastoplastic':
            s_p = C * self.EAp
            relief = self.sigma_y * C / (C + H) if self.alpha > 0.0 else 0.0
        else:
            s_p = relief = 0.0
        self.hist.append((tau, s_p, relief, tau + s_p + relief))
        return tau, s, gam, dgam

    @property
    def tau_total(self):
        """Surface tension composed of the elastic and plastic components at the last step."""
        return self.hist[-1][3] if self.hist else -self.sigma_0


def DT_capillary(gam, Delta_Sv, r):
    """Undercooling returned by brentq(f_DT, 1e-4, 1) in the original scripts (singularity of f_DT)."""
    return 2 * gam / (r * abs(Delta_Sv))
