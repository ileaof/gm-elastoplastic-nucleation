# Fig 11 do paper — a série E6.5 (800/850/900 K), rótulos em inglês, estilo do manuscrito.
# Dados: CSVs das análises + logs de thermo + dumps dos controles puros (pasta do laptop).
# Saída: Fig11_E65_series.png dentro de paper/ (repo OneDrive).
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
SAIDA = r"C:\Users\ileao\OneDrive\Documentos\gm-elastoplastic-nucleation\paper\Fig11_E65_series.png"
N, VOL = 62500, 1.034e-24
TEMPOS = ((800, "e65_analise_T800", ("ctrl_fcc_T800_hold0.dump", "ctrl_fcc_T800_late.dump")),
          (850, "e65_analise_T850", ("ctrl_fcc_T850_hold0.dump", "ctrl_fcc_T850_late.dump")),
          (900, "e65_analise", ("ctrl_fcc_T900.dump", "ctrl_fcc_T900_late.dump")))

def dump(path):
    lines = open(path).read().split("\n")
    n = int(lines[3])
    return np.array([l.split() for l in lines[9:9 + n]], dtype=float)

def thermo(path):
    rows = []
    for ln in open(path, errors="ignore"):
        p = ln.split()
        if len(p) == 6:
            try:
                rows.append([float(x) for x in p])
            except ValueError:
                pass
    return np.array(rows)

CORES = {800: "#29627e", 850: "#c87a2e", 900: "#8a2f2f"}
fig, axs = plt.subplots(1, 3, figsize=(13.2, 4.2))

print(f"{'T':>5} {'f5 gross':>9} {'baseline':>9} {'excess':>8} {'PE sheared':>11} "
      f"{'PE pure':>9} {'delta':>8} {'J/m3':>10} {'Pxy [kbar]':>10}")
for T, pasta, ctrls in TEMPOS:
    f5_p, pe_p = [], []
    for nome in ctrls:
        a = dump(os.path.join(AQUI, nome))
        m5 = a[:, 6].astype(int) == 5
        f5_p.append(m5.mean()); pe_p.append(a[:, 5].mean())
    f5_puro = np.mean(f5_p); pe_puro = np.mean(pe_p)
    c = np.genfromtxt(os.path.join(AQUI, pasta, "e65_resumo.csv"), delimiter=",", names=True)
    t, f5, pe = c["t_ps"], c["f_diamante"], c["PE_medio"]
    d = pe - pe_puro
    dens = d * N * 1.602176634e-19 / VOL
    lg = "rec_gpu.log" if T == 900 else f"rec_gpu_T{T}.log"
    lg_path = os.path.join(AQUI, "e65_analise", lg) if T == 900 else os.path.join(AQUI, lg)
    h = thermo(lg_path); h = h[h[:, 0] >= 15000]
    tv = (h[:, 0] - 15000) * 0.002
    print(f"{T:4d}K {100*f5.mean():8.1f}% {100*f5_puro:8.1f}% {100*(f5.mean()-f5_puro):7.1f}% "
          f"{pe.mean():11.4f} {pe_puro:9.4f} {d.mean():+8.4f} {dens.mean():9.2e} "
          f"{h[:, 4].mean()/1000:10.1f}")
    a0 = axs[0]
    a0.plot(t, pe, "o-", ms=3, color=CORES[T], label=f"{T} K, $\\gamma$ = 0.10 sheared")
    a0.plot((t[0], t[-1]), (pe_puro, pe_puro), "--", lw=1.8, color=CORES[T], alpha=0.65,
            label=f"{T} K, pure-FCC control")
    a1 = axs[1]
    a1.plot(t, 100 * f5, "o-", ms=3, color=CORES[T], label=f"{T} K")
    a1.fill_between(t, 100 * f5_puro, 100 * f5, color=CORES[T], alpha=0.16)
    a2 = axs[2]
    a2.plot(tv, h[:, 4] / 1000, ".", ms=2, color=CORES[T],
            label=f"{T} K ($\\sigma_{{xy}}$ = {h[:, 4].mean() / 1000:.0f} kbar)")

axs[0].set_ylabel("Potential energy per atom [eV]")
axs[0].set_title("(a) Energy during the 1 ns anneals\n(dashed = pristine-FCC control)", fontsize=9)
axs[1].set_ylabel("CNA label '5' fraction [%]")
axs[1].set_title("(b) Label '5': shaded area = excess over\nthe per-temperature pristine baseline", fontsize=9)
axs[2].set_ylabel("$\\sigma_{xy}$ [kbar]")
axs[2].set_title("(c) Residual shear stress\n(does not relax within 1 ns)", fontsize=9)
for ax in axs:
    ax.set_xlabel("Time in the 800-900 K anneal [ps]")
    ax.tick_params(labelsize=8.5)
    ax.grid(True, ls="--", color="0.85")
    ax.legend(fontsize=7.5, loc="best")
fig.tight_layout()
fig.savefig(SAIDA, dpi=200)
print("OK:", SAIDA)