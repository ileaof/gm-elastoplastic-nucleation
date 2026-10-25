"""Figures of the recrystallization section (5.7) of the elasto-plastic surface paper.

Fig7_E6_stored.png: stored energy E_st(eps) and shear stress sigma(eps) from the MD shear steps
  (Jelinek MEAM of the Al-Si-Mg-Fe alloy; Mendelev EAM pure Al as the control), with the model's
  yield strain E_A^p = sigma_y/|sigma_0| = 0.144 marked. Data: the GPU machine run (clean 0.30).
Fig8_E6_anneal.png: potential energy along the annealing steps: recovery at 800-900 K (control);
  the violent threshold of the alloy between 700 and 800 K (laptop runs, 2 048 atoms).
Fig9_E6_field.png: the stored-deformation FIELD as atom maps (energy per atom, colored):
  the deformed state (the defect network glows) and the annealed states (released or not).
Usage:  python make_e6_figures.py    (reads ../recristalizacao/ and ../recristalizacao/gpu/results_e6/)
"""
import os
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R6 = os.path.join(HERE, '..', 'recristalizacao')
GPU = os.path.join(R6, 'gpu', 'results_e6')
OUT = HERE
plt.rcParams.update({'font.size': 11, 'font.family': 'serif', 'mathtext.fontset': 'dejavuserif'})
COL = {'jelinek': '#2a78d6', 'mendelev': '#eb6834'}
TS = 0.002  # ps per step


def stored(src, pot):
    f = os.path.join(src, f'stored_{pot}.dat')
    log = os.path.join(src, f'deform_{pot}.log')
    rows = [l.split() for l in open(f, encoding='utf-8', errors='replace') if l.startswith('STORED')]
    ref = float(re.search(r'^REF\s+(-[\d.]+)', open(log, encoding='utf-8', errors='replace').read(),
                          re.M).group(1))
    eps = np.array([float(r[1]) for r in rows])
    eat = np.array([float(r[2]) for r in rows]) - ref
    sig = np.array([float(r[4]) for r in rows]) * 0.1        # bar -> MPa
    return eps, eat, sig


def thermo(log):
    raw = open(log, encoding='utf-8', errors='replace').read()
    best = None
    for m in re.finditer(r'^\s*Step\s+Temp\s+Press\s+PotEng[^\n]*\n(.*?)(?:\n\s*\n|\n\s*Loop time)',
                         raw, re.M | re.S):
        rows = [r for r in m.group(1).splitlines() if re.match(r'^\s*\d+\s', r)]
        if len(rows) > (0 if best is None else len(best)):
            best = np.array([[float(x) for x in r.split()[:5]] for r in rows])
    return None if best is None else (best[:, 0] * TS, best[:, 3])


def last_atoms(path):
    lines = [l.rstrip('\n') for l in open(path, encoding='utf-8', errors='replace')]
    i, cols, n, cur = 0, None, 0, None
    while i < len(lines):
        if lines[i].startswith('ITEM: NUMBER OF ATOMS'):
            n = int(lines[i + 1]); i += 2
        elif lines[i].startswith('ITEM: BOX BOUNDS'):
            i += 4
        elif lines[i].startswith('ITEM: ATOMS'):
            cols = {c: k for k, c in enumerate(lines[i].split()[2:])}
            i += 1
            cur = np.array([[float(x) for x in lines[i + k].split()[:len(cols)]] for k in range(n)])
            i += n
        else:
            i += 1
    return cols, cur


# Fig7 — stored energy and flow in the shear steps (the model's yield strain marked)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
for pot in ('jelinek', 'mendelev'):
    eps, eat, sig = stored(GPU, pot)
    label = {'jelinek': 'Jelinek MEAM, Al-0.8Si-0.6Mg-0.2Fe', 'mendelev': 'Mendelev EAM, pure Al (control)'}[pot]
    axs[0].plot(eps, eat, 'o-', color=COL[pot], lw=1.8, ms=4.5, label=label)
    axs[1].plot(eps, sig, 'o-', color=COL[pot], lw=1.8, ms=4.5)
axs[0].axvline(0.144, color='k', lw=1.2, ls='-.',
               label=r'Model yield strain, $E_A^p=\sigma_y/|\sigma_0| = 0.144$ (Al-Si-Mg-Fe)')
axs[0].set_xlabel(r'Shear strain $\varepsilon$')
axs[0].set_ylabel(r'Stored energy per atom, $E_{st}$ [eV]')
axs[1].set_xlabel(r'Shear strain $\varepsilon$')
axs[1].set_ylabel(r'Shear stress $\sigma_{xy}$ [MPa]')
axs[0].legend(fontsize=8.5, frameon=True, edgecolor='k')
for ax in axs:
    ax.grid(True, ls='--', color='0.8', lw=0.5)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'Fig7_E6_stored.png'), dpi=300)
plt.close(fig)

# Fig8 — annealing steps: recovery (control, 800-900 K) and the violent threshold (alloy, 700 K)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
cmap = plt.get_cmap('viridis')
TJ = {'jelinek': [500, 600, 700], 'mendelev': [500, 600, 700, 800, 850, 900]}
for i, (pot, title) in enumerate((('jelinek', 'Jelinek MEAM (alloy)'),
                                  ('mendelev', 'Mendelev EAM (control)'))):
    logs = [os.path.join(R6, f'recozer_{pot}_T{t:g}.log') for t in TJ[pot]]
    for k, log in enumerate([f for f in logs if os.path.exists(f)]):
        th = thermo(log)
        if th is None:
            continue
        t_ps, pe = th
        T = float(re.search(r'_T([0-9.]+)\.log', log).group(1))
        c = cmap(k / max(len(logs) - 1, 1))
        axs[i].plot(t_ps, pe, color=c, lw=1.0, label=f'T = {T:g} K')
    axs[i].set_title(title, fontsize=11)
    axs[i].set_xlabel('time [ps]')
    axs[i].set_ylabel(r'PE per atom [eV]')
    axs[i].legend(fontsize=7.5, frameon=True, edgecolor='k', ncol=2)
    axs[i].grid(True, ls='--', color='0.8', lw=0.5)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'Fig8_E6_anneal.png'), dpi=300)
plt.close(fig)

# Fig10 — the onset times of the catastrophic decay across the fine frontier (the stochastic map)
fig = plt.figure(figsize=(10.5, 4.4))
ax = fig.add_axes([0.09, 0.16, 0.78, 0.72])
EVENTS = [
    (750, 11, 41), (750, 22, 64), (750, 33, 17),
    (775, 11, 19), (775, 33, 66),
    (800, 33, 56),
]
CLEAN = [(775, 22), (800, 11), (800, 22)]
for T, S, t in EVENTS:
    ax.scatter(T, t, c='#eb6834', s=52, zorder=5)
    ax.annotate(f's{S}', (T, t), xytext=(4, -2), textcoords='offset points', fontsize=8, color='#eb6834')
CL = {775: ['s22'], 800: ['s11', 's22']}
for T, seeds in CL.items():
    ax.scatter(T, 100, c='#2a78d6', s=52, marker='s', zorder=5)
    ax.annotate(' e '.join(seeds) + ' limpo', (T, 100), xytext=(6, -11), textcoords='offset points',
                fontsize=8, color='#2a78d6')
ax.axhline(100, color='k', lw=1.2)
ax.set_xlabel('temperatura do recozimento [K]')
ax.set_ylabel('tempo do início do decaimento [ps]')
ax.set_xticks([750, 775, 800])
ax.grid(True, ls='--', color='0.85')
fig.text(0.5, 0.945, 'The onset of the catastrophic decay at 750-800 K: 6 of 9 anneals decayed within 17-66 ps',
         ha='center', fontsize=11)
fig.savefig(os.path.join(OUT, 'Fig10_E6_onset.png'), dpi=300)
print('Fig10 ok')

# Fig9 — the stored-deformation field: atom maps (energy per atom), deformed x annealed.
# The deformed state is the FIRST frame of the anneal dumps (the anneal restarts from the
# sheared configuration); the annealed state is its last frame.
fig = plt.figure(figsize=(13.5, 5.0))
axs = [fig.add_axes([0.03 + 0.245 * k, 0.22, 0.225, 0.52]) for k in range(4)]
panels = [
    ('recozido_mendelev_T500.dump', 'first', 'Mendelev, deformed ($\\varepsilon$ = 0.14)'),
    ('recozido_mendelev_T500.dump', 'last', 'Mendelev, annealed 500 K (still stored)'),
    ('recozido_mendelev_T800.dump', 'last', 'Mendelev, annealed 800 K (released)'),
    ('recozido_jelinek_T700.dump', 'first', 'Jelinek (alloy), deformed ($\\varepsilon$ = 0.14)'),
]


def atom_frame(path, which):
    lines = [l.rstrip('\n') for l in open(path, encoding='utf-8', errors='replace')]
    frames = []
    i, cols, n = 0, None, 0
    while i < len(lines):
        if lines[i].startswith('ITEM: NUMBER OF ATOMS'):
            n = int(lines[i + 1]); i += 2
        elif lines[i].startswith('ITEM: BOX BOUNDS'):
            i += 4
        elif lines[i].startswith('ITEM: ATOMS'):
            cols = {c: k for k, c in enumerate(lines[i].split()[2:])}
            i += 1
            frames.append(np.array([[float(x) for x in lines[i + k].split()[:len(cols)]]
                                    for k in range(n)]))
            i += n
        else:
            i += 1
    if not frames or cols is None:
        return (cols, None)[0], None
    return cols, (frames[0] if which == 'first' else frames[-1])


imgs = []
for ax, (fn, which, title) in zip(axs, panels):
    src = os.path.join(R6, fn)
    cols, cur = atom_frame(src, which) if os.path.exists(src) else (None, None)
    if cur is None:
        ax.axis('off')
        continue
    im = ax.scatter(cur[:, cols['x']], cur[:, cols['z']], s=2.0, c=cur[:, cols['c_esp']],
                    cmap='inferno', vmin=-3.42, vmax=-3.20, lw=0)
    imgs.append(im)
    ax.set_title(title, fontsize=10)
    ax.set_xticks([]), ax.set_yticks([])
cb = fig.colorbar(imgs[0], cax=fig.add_axes([0.32, 0.085, 0.36, 0.030]), orientation='horizontal')
cb.set_label(r'Energy per atom [eV]  —  bright = stored in the defect network', fontsize=10)
cb.ax.hlines(0.5, -3.42, -3.20, color='none')
cb.ax.tick_params(labelsize=9)
fig.text(0.5, 0.955, 'The stored plastic deformation as an energy field of the deformed solid, and its release by annealing',
         ha='center', fontsize=11)
fig.savefig(os.path.join(OUT, 'Fig9_E6_field.png'), dpi=300)
plt.close(fig)
print('Figs 7-9 ok')