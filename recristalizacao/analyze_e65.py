# E6.5 — análise do recozimento a 900 K da célula grande do controle Mendelev
# (62 500 átomos, cisalhamento + 1 ns, dumps a cada 25 000 passos = 50 ps).
# Por frame: frações CNA (LAMMPS: 0 other, 1 fcc, 2 hcp, 3 bcc, 4 ico, 5 d Diamante),
# energia potencial média global e por grupo; CSV + figuras.
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

T = sys.argv[1] if len(sys.argv) > 1 else "900"
DUMP = os.path.expanduser(f"~/e6_rec/reciz_mendelev_T{T}.dump")
OUTDIR = os.path.expanduser(f"~/e6_rec/e65_analise_T{T}")
STEP_ANNEAL_INICIO = 15000   # 5000 (300 K) + 5000 (cisalh.) + 5000 (rampa)
DT = 0.002                   # ps
NOMES = {0: "other", 1: "fcc", 2: "hcp", 3: "bcc", 4: "ico", 5: "diamante"}

os.makedirs(OUTDIR, exist_ok=True)

steps, tempos, fracoes, pes = [], [], [], []   # fracoes: [N,6]; pes: [N,3]
with open(DUMP, "r") as f:
    lines = f.read().split("\n")
i, n_frames = 0, 0
while i < len(lines):
    if lines[i].startswith("ITEM: TIMESTEP"):
        step = int(lines[i + 1])
        natoms = int(lines[i + 3])
        n_linhas = i + 9 + natoms
        arr = np.array([l.split() for l in lines[i + 9:n_linhas] if l.strip()], dtype=float)
        cna = arr[:, 6].astype(int)
        pe = arr[:, 5]
        f_cna = np.array([(cna == k).mean() for k in range(6)])
        pes.append([pe.mean(), pe[cna == 1].mean() if (cna == 1).any() else np.nan,
                    pe[cna != 1].mean() if (cna != 1).any() else np.nan])
        fracoes.append(f_cna)
        steps.append(step)
        tempos.append((step - STEP_ANNEAL_INICIO) * DT)
        n_atoms_frame = natoms
        n_frames += 1
        i = n_linhas
    else:
        i += 1

fracoes = np.array(fracoes)
pes = np.array(pes)
tempos = np.array(tempos)

print(f"frames: {n_frames}, átomos/frame: {n_atoms_frame}")
print(f"{'t_ps':>7} {'fcc':>6} {'hcp':>6} {'other+':>7} {'diam':>6} {'PE_med':>8} {'PE_fcc':>8} {'PE_ndef':>8}")
for j, t in enumerate(tempos):
    print(f"{t:7.0f} {fracoes[j,1]:6.3f} {fracoes[j,2]:6.3f} {fracoes[j,0]+fracoes[j,3]+fracoes[j,4]:7.3f} "
          f"{fracoes[j,5]:6.3f} {pes[j,0]:8.3f} {pes[j,1]:8.3f} {pes[j,2]:8.3f}")

cab = ("t_ps,step,f_fcc,f_hcp,f_bcc,f_ico,f_diamante,f_other,PE_medio,PE_fcc,PE_nao_fcc")
with open(os.path.join(OUTDIR, "e65_resumo.csv"), "w") as g:
    g.write(cab + "\n")
    for j, t in enumerate(tempos):
        g.write(f"{t:.0f},{steps[j]},{fracoes[j,1]:.6f},{fracoes[j,2]:.6f},{fracoes[j,3]:.6f},"
                f"{fracoes[j,4]:.6f},{fracoes[j,5]:.6f},{fracoes[j,0]:.6f},{pes[j,0]:.6f},"
                f"{pes[j,1]:.6f},{pes[j,2]:.6f}\n")

# Figura: 2 painéis
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
a = ax[0]
a.plot(tempos, 100 * fracoes[:, 1], "o-", color="#2f6f6f", label="FCC (CNA 1)")
a.plot(tempos, 100 * fracoes[:, 2], "s-", color="#c87a2e", label="HCP (CNA 2)")
resto = 100 * (fracoes[:, 0] + fracoes[:, 3] + fracoes[:, 4])
a.plot(tempos, resto, "^-", color="#5a4a82", label="defeitos (other+bcc+ico)")
a.plot(tempos, 100 * fracoes[:, 5], "d-", color="#8a8a8a", label="diamond (CNA 5)", ms=4)
a.set_xlabel("tempo no recozimento 900 K (ps)")
a.set_ylabel("fração de átomos (%)")
a.set_title(f"E6.5 — Mendelev 62 500 átomos, {T} K, 1 ns\nfrações CNA por frame", fontsize=10)
a.legend(fontsize=8)
a.grid(alpha=0.3)
b = ax[1]
b.plot(tempos, pes[:, 0], "o-", color="#1c3f5a", label="todos os átomos")
b.plot(tempos, pes[:, 1], "s-", color="#2f6f6f", label="átomos FCC")
b.plot(tempos, pes[:, 2], "^-", color="#c87a2e", label="átomos não-FCC")
b.set_xlabel("tempo no recozimento 900 K (ps)")
b.set_ylabel("E_pot por átomo (eV)")
b.set_title("E6.5 — energia potencial média por grupo CNA", fontsize=10)
b.legend(fontsize=8)
b.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, f"fig_E65_T{T}.png"), dpi=150)
print("CSV e figura em:", OUTDIR)