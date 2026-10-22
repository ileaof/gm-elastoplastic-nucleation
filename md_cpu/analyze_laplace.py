"""
E5 - interface stress of a critical seed from its Laplace pressure.

    f = dP * r / 2,   dP = P_solid(core sphere) - P_liquid(far region)

compared with the interfacial energy of the same seed, gamma = |dS_V| * Gamma (JMEP Eq. 5,
Table 4), so that the Shuttleworth excess  f - gamma  (Gurtin-Murdoch excess surface stress)
can be read for each seed size.

Usage:  python analyze_laplace.py laplace_r20.0_T809.3_n24.dat --r0 20 --gamma-ref 1.217e-7
"""
import argparse
import numpy as np

dSv_MD = 8.878e5            # |dS_V| Borovikov EAM [J m^-3 K^-1] (E2_results.md)
BAR = 1.0e5                 # Pa


def load(fname):
    a = np.loadtxt(fname, comments='#')
    return {'step': a[:, 0], 'Ns': a[:, 1], 'Pin': a[:, 2], 'Pout': a[:, 3], 'dP': a[:, 4],
            'T': a[:, 5], 'P': a[:, 6]}


def stationary_window(Ns, tol=0.15, skip=4):
    """Samples after `skip` (equilibration of the released seed) whose solid count stays within
    +/- tol of its value at the start of the window."""
    N_ref = np.median(Ns[skip:skip + 4])
    ok = np.abs(Ns - N_ref) <= tol * N_ref
    ok[:skip] = False
    # keep the first contiguous block only (seed growing or melting afterwards is excluded)
    idx = np.where(ok)[0]
    if idx.size == 0:
        return ok
    end = idx[0]
    while end + 1 < len(ok) and ok[end + 1]:
        end += 1
    w = np.zeros_like(ok)
    w[idx[0]:end + 1] = True
    return w


def block_error(x, nblocks=5):
    n = len(x) // nblocks
    if n < 1:
        return float('nan')
    means = [x[k * n:(k + 1) * n].mean() for k in range(nblocks)]
    return np.std(means, ddof=1) / np.sqrt(nblocks)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('file')
    ap.add_argument('--r0', type=float, required=True, help='seed radius at the start of production [A]')
    ap.add_argument('--gamma-ref', type=float, default=None, help='Gamma of this seed (JMEP Table 4) [m K]')
    ap.add_argument('--dt', type=float, default=0.002)
    args = ap.parse_args()

    d = load(args.file)
    w = stationary_window(d['Ns'])
    t = d['step'] * args.dt
    print(f'file: {args.file}')
    print(f'samples: {len(t)}, stationary window: {w.sum()} samples, '
          f't = {t[w][0]:.1f}-{t[w][-1]:.1f} ps' if w.sum() else 'no stationary window')
    print(f'N_solid: start {d["Ns"][0]:.0f}, window mean {d["Ns"][w].mean():.0f}, end {d["Ns"][-1]:.0f}')
    print(f'T = {d["T"][w].mean():.1f} K,  box P = {d["P"][w].mean():.0f} bar')
    if w.sum() < 10:
        print('WARNING: window too short for a reliable pressure average')
    dP = d['dP'][w].mean()
    dP_err = block_error(d['dP'][w])
    r = args.r0 * 1e-10
    f = dP * BAR * r / 2
    f_err = dP_err * BAR * r / 2
    print(f'P_in = {d["Pin"][w].mean():.0f} bar, P_out = {d["Pout"][w].mean():.0f} bar, '
          f'dP = {dP:.0f} +/- {dP_err:.0f} bar')
    print(f'f = dP r / 2 = {f:.4f} +/- {f_err:.4f} N/m   (r = {args.r0:.1f} A)')
    if args.gamma_ref:
        gam = dSv_MD * args.gamma_ref
        print(f'gamma = |dS_V| Gamma = {gam:.4f} J/m^2  (Gamma = {args.gamma_ref:.4e} m K)')
        print(f'Shuttleworth excess  f - gamma = {f - gam:+.4f} +/- {f_err:.4f} N/m')


if __name__ == '__main__':
    main()
