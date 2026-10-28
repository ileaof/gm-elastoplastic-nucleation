# E6.5 série 800/850/900 K — figura final no laptop.
# Baselines: dos dumps de controle baixados (hold0 + late por T);
# E6.5: dos CSVs das análises; Pxy: dos logs de thermo.
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
N, VOL = 62500, 1.034e-24
TS = ((800, "e65_analise_T800", None), (850, "e65_analise_T850", None),
      (900, "e65_analise", "e65_analise"))   # (T, pasta do csv da corrida, pasta do log)

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

CORES = {800: "#2f6f8f", 850: "#c87a2e", 900: "#8a2f2f"}
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6))
print(f"{'T':>5} {'f5_bruto':>9} {'baseline':>9} {'excesso':>8} {'PE_cisi':>9} {'PE_puro':>9} {'delta':>8} {'J/m3':>10}")

for T, pasta_csv, pasta_log in TS:
    pasta_log = pasta_log or pasta_csv
    # baseline puro: média hold0+late
    if T == 900:
        bs = ["ctrl_fcc_T900.dump", "ctrl_fcc_T900_late.dump"]
    else:
        bs = [f"ctrl_fcc_T{T}_hold0.dump", f"ctrl_fcc_T{T}_late.dump"]
    f5_p, pe_p = [], []
    for b in bs:
        a = dump(os.path.join(AQUI, b))
        m5 = a[:, 6].astype(int) == 5
        f5_p.append(m5.mean()); pe_p.append(a[:, 5].mean())
    f5_puro = float(np.mean(f5_p)); pe_puro = float(np.mean(pe_p))
    c = np.genfromtxt(os.path.join(AQUI, pasta_csv, "e65_resumo.csv"), delimiter=",", names=True)
    t = c["t_ps"]; f5 = c["f_diamante"]; pe = c["PE_medio"]
    d = pe - pe_puro
    dens = d * N * 1.602176634e-19 / VOL
    print(f"{T:4d}K {100 * f5.mean():8.1f}% {100 * f5_puro:8.1f}% {100 * (f5.mean() - f5_puro):7.1f}% "
          f"{pe.mean():9.4f} {pe_puro:9.4f} {d.mean():+8.4f} {dens.mean():9.2e}")
    # log (Pxy): 900 = e65_analise/rec_gpu.log; 800/850 = e65_analise_T*/rec_gpu_T{T}.log
    lg = "rec_gpu_T800.log" if T == 800 else ("rec_gpu_T850.log" if T == 850 else "rec_gpu.log")
    pasta_lg = os.path.join(AQUI, "e65_analise") if T == 900 else AQUI
    h = thermo(os.path.join(pasta_lg, lg))
    h = h[h[:, 0] >= 15000]
    tv = (h[:, 0] - 15000) * 0.002
    a0 = ax[0]; a0.plot(t, pe, "o-", ms=3, color=CORES[T], label=f"{T} K (ΔPE +{d.mean():.4f})")
    a0.plot((tv[0], tv[-1]), (pe_puro, pe_puro), "--", lw=2, color=CORES[T], alpha=0.6,
            label=f"puro {T} K")
    a1 = ax[1]
    a1.plot(t, 100 * f5, "o-", ms=3, color=CORES[T], label=f"{T} K")
    a1.fill_between(t, 100 * f5_puro, 100 * f5, color=CORES[T], alpha=0.15)
    a2 = ax[2]; a2.plot(tv, h[:, 4] / 1000, ".", ms=2, color=CORES[T],
                        label=f"{T} K (−{abs(h[:, 4].mean()) / 1000:.1f} kbar)")
ax[0].set_ylabel("E_pot por átomo (eV)")
ax[0].set_title("PE/átomo no hold (traço = puro)", fontsize=10)
ax[1].set_ylabel("fração rótulo '5' (%)")
ax[1].set_title("rótulo CNA '5': a área = excesso sobre\no baseline puro (17,7/22,8/28,6%)", fontsize=10)
ax[2].set_ylabel("Pxy (kbar)")
ax[2].set_title("tensão shear residual (não relaxa)", fontsize=10)
for a in ax:
    a.set_xlabel("tempo no recozimento (ps)")
    a.grid(alpha=0.3)
    a.legend(fontsize=7.5)
fig.suptitle("E6.5 — Mendelev 62 500 átomos, γ = 0,10 sem descarregamento, NVT 1 ns: "
             "o estado cisalhado fica meta-estável em 800–900 K (ΔPE ≈ +0,026 eV/átomo, "
             "≈ 2,5×10⁸ J/m³, sem liberação)", fontsize=10.5, y=1.03)
fig.tight_layout()
fig.savefig(os.path.join(AQUI, "e65_analise", "fig_E65_serie.png"), dpi=150, bbox_inches="tight")
print("figura: e65_analise/fig_E65_serie.png")