# E6.6 — figura-síntese local (2×3): PE(t) dos holds contra as réguas (linha de cima) e
# excesso do rótulo '5' sobre a régua correta por família (linha de baixo), em 800/850/900 K.
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
N, VOLAT = 62500, 1.655e-29   # volume por átomo [m3]
C = {"800": "#29627e", "850": "#c87a2e", "900": "#8a2f2f"}

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

def hold(path, passo_ini):
    h = thermo(path)
    return h[h[:, 0] >= passo_ini]

# réguas
reg = {}
for ln in open(os.path.join(AQUI, "e66_analise", "reguas.csv")):
    if ln.startswith("#") or ln.startswith("tipo") or not ln.strip():
        continue
    w = ln.split(",")
    reg[(w[0], int(w[1]))] = (float(w[2]), float(w[3]))   # (f5, PE dump)

def t_ps(h, ini):
    return (h[:, 0] - ini) * 0.002

INIS = {"largo0.20": 20000, "largo0.30": 25000}
fig, axs = plt.subplots(2, 3, figsize=(16.5, 9.0), sharex="col")

for j, T in enumerate((800, 850, 900)):
    a = axs[0, j]; b = axs[1, j]
    # ---- famílias NPT (descarregadas): γ=0.10 (média das 2 sementes) e γ=0.30
    # nomes reais: desca_npt_T{T}_s{seed}.log / desca_npt30_T{T}.log
    for gam, arqs, cor, lb in (
            ("0.10", [f"desca_npt_T{T}_s101.log", f"desca_npt_T{T}_s202.log"], "#2e5340",
             f"$\\gamma$=0.10 unloaded + NPT (2 seeds)"),
            ("0.30", [f"desca_npt30_T{T}.log"], "#7a4023",
             "$\\gamma$=0.30 unloaded + NPT")):
        pess = []
        tt = None
        for nome in arqs:
            cr = [ln for ln in open(os.path.join(AQUI, nome), errors="replace")
                  if ln.startswith("DESCA_CROSS")][0]
            pc = int(float(re.search(r"step=(\d+)", cr).group(1)))
            h = hold(os.path.join(AQUI, nome), pc + 7000)
            if tt is None:
                tt = t_ps(h, pc + 7000)
            pess.append(h[:, 3] / N)
        pess = np.mean(pess, axis=0)
        a.plot(tt, pess, "-", lw=1.4, color=cor, label=lb)
        # excesso do rótulo 5 (frames)
        fcsv = (AQUI + f"/e66_analise/frames_descarnpt_T{T}_g{gam}_s101.csv"
                if gam == "0.30" else
                "/".join([AQUI, "e66_analise", f"frames_descarnpt_T{T}_g{gam}_s202.csv"]))
        if os.path.exists(fcsv):
            fr = np.genfromtxt(fcsv, delimiter=",", names=True)
            base = reg[("NPT", T)][0]
            b.plot(fr["t_ps"], 100 * (fr["f5"] - base), "o-", ms=3, color=cor)
        a.axhline(reg[("NPT", T)][1], ls="--", lw=1.5, color=cor, alpha=0.5)
        a.axhline(reg[("NVT", T)][1], ls=":", lw=1.5, color="#555", alpha=0.7)
        b.axhline(0.0, ls="--", lw=1.2, color="#555", alpha=0.6)
    # ---- famílias NVT (carregadas): largas 0.20/0.30 + régua NVT (puro)
    for gam in ("0.20", "0.30"):
        h = hold(os.path.join(AQUI, f"rec_largo{int(round(float(gam) * 100))}_T{T}.log"),
                 INIS["largo" + gam])
        cor = "#4a7fa5" if gam == "0.20" else "#d95f02"
        a.plot(t_ps(h, INIS["largo" + gam]), h[:, 3] / N, "-", lw=1.1, color=cor,
               alpha=0.85, label=f"$\\gamma$={float(gam):.1f} loaded + NVT (Pxy ${h[:, 4].mean() / 1000:.0f}$ kbar)")
        fc = os.path.join(AQUI, "e66_analise", f"frames_largo_T{T}_g{gam}.csv")
        if os.path.exists(fc):
            fr = np.genfromtxt(fc, delimiter=",", names=True)
            base = reg[("NVT", T)][0]
            b.plot(fr["t_ps"], 100 * (fr["f5"] - base), "s-", ms=3, color=cor)
    reg_ = ax = None
    a.axhline(reg[("NVT", T)][1], ls=":", lw=1.5, color="#444", label="pure (NVT) / pure (NPT) rulers")
    a.set_title(f"{T} K", fontsize=11)
    b.set_title(f"label '5' excess over ruler (%)", fontsize=9)
    b.set_ylim(-3, 32)
    if j == 0:
        a.set_ylabel("Potential energy per atom [eV]")
        b.set_ylabel("CNA '5' excess [%]")
for x in axs[:, 0]:
    x.set_xlabel("Time in the anneal [ps]")
for a in axs[0]:
    a.legend(fontsize=7, loc="lower right")
fig.suptitle("E6.6 — 62 500-atom Mendelev cells: unloaded states annealed NPT (P = 0) vs the loaded "
             "E6.5 design (NVT); two rulers = the pristine control of each barostat; "
             "zero lost atoms anywhere", fontsize=10.5, y=0.995)
fig.tight_layout()
fig.savefig(os.path.join(AQUI, "e66_analise", "fig_E66_serie.png"), dpi=150, bbox_inches="tight")

# ---- desconto-síntese do descarregamento (do descarrega.csv)
print("\n--- descarrega.csv (resumo) ---")
d = np.genfromtxt(os.path.join(AQUI, "e66_analise", "descarrega.csv"), delimiter=",", names=True,
                  dtype=None, encoding="utf-8")
print(f"{'run':24s} {'gam':>4} {'e_cruz':>7} {'resto eV/at':>11} {'J/m2':>6} {'dPE hold':>9}")
for r in d:
    print(f"{str(r['tag']):24s} {float(r['gam']):4.2f} {float(r['eps_cross']):7.3f} "
          f"{float(r['resto_eV_at']):11.5f} {float(r['resto_Jm2']):6.2f} "
          f"{float(r['dPE_fim']):+9.5f}")
print("\nfigura:", os.path.join(AQUI, "e66_analise", "fig_E66_serie.png"))