"""
E5 - planar solid-liquid interface stress from the z-profile of in.planar_stress.

D(z) = P_zz(z) - (P_xx(z)+P_yy(z))/2, with P_k(z) = -n(z) <s_kk>(z) (per-atom virial, bar).
Excess over the bulk phases with equimolar Gibbs dividing surfaces (two interfaces):
    f = (1/2) [ int D dz - L_s D_s - L_l D_l ],   L_s = (N - n_l L) / (n_s - n_l),  L_l = L - L_s
where D_s, D_l, n_s, n_l are averages over bulk bins (>= `margin` A from the interfaces).
"""
import numpy as np

BAR_A_TO_N_M = 1.0e-5


def read_blocks(fname):
    blocks, cur = [], None
    for line in open(fname):
        if line.startswith('#'):
            continue
        p = line.split()
        if len(p) == 3:
            if cur:
                blocks.append(np.array(cur))
            cur = []
        else:
            cur.append([float(x) for x in p[1:]])
    if cur:
        blocks.append(np.array(cur))
    return blocks          # columns: z, Ncount, n, s_xx, s_yy, s_zz


def smooth(x, w=5):
    k = np.ones(w) / w
    return np.convolve(np.concatenate([x[-w:], x, x[:w]]), k, mode='same')[w:-w]


def interface_stress(b, margin=6.0, trim=2, w=9):
    b = b[trim:len(b) - trim]          # edge bins are partially outside the (NPH) box
    z, n = b[:, 0], b[:, 2]
    dz = z[1] - z[0]
    L = dz * len(z)
    P = -b[:, 3:6] * n[:, None]
    D = P[:, 2] - 0.5 * (P[:, 0] + P[:, 1])
    ns = smooth(n, w)
    thr = 0.5 * (ns.max() + ns.min())
    sol = ns > thr
    # bins far from any solid/liquid switch
    switches = np.where(sol != np.roll(sol, 1))[0]
    zc = z[switches]
    dist = np.min(np.abs(((z[:, None] - zc[None, :]) + L / 2) % L - L / 2), axis=1)
    bulk = dist >= margin
    n_s, n_l = n[sol & bulk].mean(), n[~sol & bulk].mean()
    D_s, D_l = D[sol & bulk].mean(), D[~sol & bulk].mean()
    N = (n * dz).sum()
    L_s = (N - n_l * L) / (n_s - n_l)
    f = 0.5 * ((D * dz).sum() - L_s * D_s - (L - L_s) * D_l) * BAR_A_TO_N_M
    return f, dict(n_s=n_s, n_l=n_l, D_s=D_s, D_l=D_l, L_s=L_s, L=L, ninterf=len(switches))


if __name__ == '__main__':
    import sys
    fname = sys.argv[1] if len(sys.argv) > 1 else 'planar_profile.dat'
    fs = []
    for k, b in enumerate(read_blocks(fname)):
        f, info = interface_stress(b)
        fs.append(f)
        print(f'block {k+1}: f = {f:+.4f} N/m | n_s {info["n_s"]:.4f}, n_l {info["n_l"]:.4f} A^-3, '
              f'D_s {info["D_s"]:+.0f}, D_l {info["D_l"]:+.0f} bar, L_s/L {info["L_s"]/info["L"]:.3f}, '
              f'switches {info["ninterf"]}')
    fs = np.array(fs[1:])      # first block = equilibration
    print(f'f_planar = {fs.mean():+.4f} +/- {fs.std(ddof=1)/np.sqrt(len(fs)):.4f} N/m  (blocks 2-{len(fs)+1})')
