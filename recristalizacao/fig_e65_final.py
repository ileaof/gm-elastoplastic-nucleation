# E6.5 — figura final da análise (no laptop): PE, fração CNA "5" com baseline do
# controle puro e Pxy residual. Entradas: e65_analise/e65_resumo.csv e rec_gpu.log.
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
DT = 0.002  # ps
NAT = 62500

csv = np.genfromtxt(os.path.join(AQUI, "e65_analise", "e65_resumo.csv"), delimiter=",", names=True)
t = csv["t_ps"]
f5 = csv["f_diamante"]          # rótulo "5" do cna/atom (bruto)
pe = csv["PE_medio"]

# linha do thermo (rec_gpu.log): step T press PE Pxy volume
rows = []
with open(os.path.join(AQUI, "e65_analise", "rec_gpu.log"), errors="ignore") as g:
    for ln in g:
        p = ln.split()
        if len(p) == 6:
            try:
                rows.append([float(x) for x in p])
            except ValueError:
                pass
rows = np.array(rows)
hold = rows[rows[:, 0] >= 15000]
t_th = (hold[:, 0] - 15000) * DT

# baselines do controle puro (ctrl_fcc.log: PE média/átomo no hold)
pe_puro_900 = -3.2867
f5_puro_900 = 0.286          # fração "5" no FCC puro a 900 K (média early/late ~0.2859)
excesso = f5 - f5_puro_900
d_pe = pe - pe_puro_900      # eV/átomo
# densidade de energia armazenada: célula 1 034 131 Å3 = 1.034e-24 m3
vol_m3 = 1034131.1 * 1e-30
dens = d_pe * NAT * 1.602176634e-19 / vol_m3   # J/m3 (do PE médio por átomo)

fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6))
a = ax[0]
a.plot(t, pe, "o-", color="#1c3f5a", ms=4, label="E6.5 cisalhado (900 K)")
a.axhline(pe_puro_900, ls="--", color="#8a2f2f", label="FCC puro 900 K (controle)")
a.set_xlabel("tempo no recozimento (ps)")
a.set_ylabel("E_pot por átomo (eV)")
a.set_title(f"E6.5 — PE/átomo: Δ médio +{d_pe.mean():.4f} eV/átomo\n"
            f"(≈ +{dens.mean():.2e} J/m³ acima do puro)", fontsize=10)
a.legend(fontsize=8)
a.grid(alpha=0.3)

b = ax[1]
b.plot(t, 100 * f5, "o-", color="#5a4a82", ms=4, label=r"rótulo CNA '5' bruto")
b.axhline(100 * f5_puro_900, ls="--", color="#8a2f2f",
          label="baseline FCC puro a 900 K (28,6%)")
b.fill_between(t, 100 * f5_puro_900, 100 * f5, color="#5a4a82", alpha=0.18,
               label="excesso = sítios defeituosos reais")
b.set_xlabel("tempo no recozimento (ps)")
b.set_ylabel("fração de átomos (%)")
b.set_title(f"E6.5 — rótulo '5' do cna/atom: bruto − baseline\n"
            f"(excesso médio {100 * excesso.mean():.1f}%; artefato térmico a 900 K)", fontsize=10)
b.legend(fontsize=8)
b.grid(alpha=0.3)

c = ax[2]
c.plot(t_th, hold[:, 4] / 1000, ".", color="#c87a2e", ms=3)
c.set_xlabel("tempo no recozimento (ps)")
c.set_ylabel("Pxy (kbar)")
c.set_title(f"E6.5 — tensão de cisalhamento residual: Pxy médio "
            f"{hold[:, 4].mean() / 1000:.1f} kbar, não relaxa em 1 ns", fontsize=10)
c.grid(alpha=0.3)

fig.suptitle("E6.5 — Mendelev 62 500 átomos, cisalhamento γ = 0,10 (tilt), rampa blindada, "
             "recoximento NVT 900 K × 1 ns: estado meta-estável, sem recristalização", fontsize=11, y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(AQUI, "fig_E65_resumo.png"), dpi=150, bbox_inches="tight")
print("OK — fig_E65_resumo.png em", AQUI)
print(f"PE medio E6.5 {pe.mean():.4f}; puro {pe_puro_900}; delta {d_pe.mean():+.4f} eV/atom; "
      f"densidade {dens.mean():.3e} J/m3; excesso '5' {100 * excesso.mean():.1f}%")