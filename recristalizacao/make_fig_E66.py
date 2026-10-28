# Fig 12 do paper (E6.6), versão limpa:
# (a) σxy(ε) em 300 K: o ramo elástico compartilhado (verde, ida+volta em ε=0,00) seguido pelo
#     ramo do γ=0.30 (marrom: carga até −73 kbar/0.18, deslizamento, descarregamento cruzando em 0.252);
# (b-d) o excesso de PE sobre a régua do barostato próprio, bruta fina + média móvel, % liberado em caixa.
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
SAIDA = r"C:\Users\ileao\OneDrive\Documentos\gm-elastoplastic-nucleation\paper\Fig12_E66_unload_release.png"
N, DT = 62500, 0.002

def thermo(caminho):
    rows = []
    for ln in open(caminho, errors="ignore"):
        p = ln.split()
        if len(p) == 6:
            try:
                rows.append([float(x) for x in p])
            except ValueError:
                pass
    return np.array(rows)

def leg(caminho, ini):
    h = thermo(caminho)
    return h[h[:, 0] >= ini]

def media_movel(y, k):   # (unused)
    return suaviza(y, k)

def suaviza(y, k):
    """Média móvel com padding de borda: mesma dimensão de y, sem cap de degrau."""
    y = np.asarray(y, float)
    z = np.pad(y, (k // 2, k - 1 - k // 2), mode="edge")
    return np.convolve(z, np.ones(k) / k, "valid")

def cross(pasta, nome):
    cr = [ln for ln in open(os.path.join(pasta, nome), errors="replace") if ln.startswith("DESCA_CROSS")][0]
    return int(float(re.search(r"step=(\d+)", cr).group(1)))

# réguas
r_npt = {}
for T in (800, 850, 900):
    h = leg(os.path.join(AQUI, f"ctrl_npt_T{T}.log"), 10002)
    r_npt[T] = h[:, 3].mean() / N
r_nvt = {}
for ln in open(r"C:\Users\ileao\OneDrive\Documentos\gm-elastoplastic-nucleation"
               r"\recristalizacao\e65_serie26\baseline_puro.dat"):
    if not (ln.startswith("#") or not ln.strip()):
        w = ln.split()
        r_nvt[int(w[0])] = float(w[2])

fig = plt.figure(figsize=(14.5, 6.0))
gs = fig.add_gridspec(1, 4, width_ratios=(1.15, 1, 1, 1), wspace=0.40)
axa = fig.add_subplot(gs[0, 0])

VERDE, MARROM = "#2e6b3f", "#a1501f"

# ---- (a)
h10 = thermo(os.path.join(AQUI, "desca_npt_T800_s101.log"))
ida10 = h10[(h10[:, 0] >= 5000) & (h10[:, 0] <= 10001)]
cr10 = cross(AQUI, "desca_npt_T800_s101.log")
vol10 = h10[(h10[:, 0] >= 10000) & (h10[:, 0] <= cr10 + 1)]
axa.plot(2e-5 * (ida10[:, 0] - 5000), ida10[:, 4] / 1000, "-", lw=1.8, color=VERDE)
axa.plot(0.10 - 2e-5 * (vol10[:, 0] - 10000), vol10[:, 4] / 1000, "--", lw=1.5, color=VERDE)
axa.plot(0.10 - 2e-5 * (cr10 - 10000), 47.07 / 1000, "o", ms=8, mfc="none", color=VERDE)
h30 = thermo(os.path.join(AQUI, "desca_npt30_T800.log"))
ida30 = h30[(h30[:, 0] >= 10000) & (h30[:, 0] <= 20001)]
cr30 = cross(AQUI, "desca_npt30_T800.log")
vol30 = h30[(h30[:, 0] >= 20000) & (h30[:, 0] <= cr30 + 1)]
axa.plot(2e-5 * (ida30[:, 0] - 5000), ida30[:, 4] / 1000, "-", lw=1.8, color=MARROM)
axa.plot(0.30 - 2e-5 * (vol30[:, 0] - 20000), vol30[:, 4] / 1000, "--", lw=1.5, color=MARROM)
axa.annotate("$\\gamma$=0.10: the unload retraces\nthe same path to $\\varepsilon$=0.00\n(remainder $\\approx$ 0)",
             (0.075, -13), fontsize=8, color=VERDE, ha="center")
axa.annotate("glide", (0.208, -33), fontsize=9, color=MARROM, style="italic")
axa.plot([0.252], [0.237], "s", ms=8, mfc="none", color=MARROM)
axa.annotate("crosses at $\\varepsilon$=0.252", (0.252, 0.237), fontsize=8, color=MARROM,
             ha="center", xytext=(0.205, 4.5), textcoords="data",
             arrowprops=dict(arrowstyle="->", color=MARROM, lw=0.8))
axa.set_xlabel("shear strain $\\varepsilon$")
axa.set_ylabel("$\\sigma_{xy}$ [kbar]")
axa.set_title("(a) shear and unloading at 300 K", fontsize=10)
axa.set_xlim(-0.02, 0.33)
axa.set_ylim(-78, 8)
axa.legend(handles=[plt.Line2D([], [], color=VERDE, lw=1.8, label="elastic load/unload ($\\gamma$=0.10)"),
                    plt.Line2D([], [], color=MARROM, lw=1.8, label="load through the glide ($\\gamma$=0.30)")],
           fontsize=7.5, loc="lower left")
axa.grid(True, ls="--", color="0.85")

# ---- (b-d)
FAMI = (("$\\gamma$=0.30 unloaded + NPT", "desca_npt30_T{T}.log", "NPT", "#5a2d1a", "o", 29400),
        ("$\\gamma$=0.20 loaded + NVT", "rec_largo20_T{T}.log", "NVT", "#3a6ea5", "s", 20000),
        ("$\\gamma$=0.30 loaded + NVT", "rec_largo30_T{T}.log", "NVT", "#d95f02", "^", 25000))
for k, T in enumerate((800, 850, 900)):
    ax = fig.add_subplot(gs[0, k + 1])
    lib_txt = []
    for lb, arq, chave, cor, mk, base in FAMI:
        h = leg(os.path.join(AQUI, arq.format(T=T)), base)
        r = r_npt if chave == "NPT" else r_nvt
        exc = (h[:, 3] / N - r[T]) * 1000.0
        tt = (h[:, 0] - base) * DT
        ax.plot(tt, exc, "-", lw=0.6, color=cor, alpha=0.30)
        ax.plot(tt, suaviza(exc, 41), "-", lw=1.7, color=cor)
        ini, fim = exc[0], exc[-1]
        lib = 100 * (1 - fim / ini) if ini > 5.0 else float("nan")
        lib_txt.append((cor, f"{lb}: {lib:.0f} %"))
        print(f"(bd) T={T} {lb}: exc {ini:+.1f} -> {fim:+.1f} meV/at; libera {lib:.0f}%")
    ax.axhline(0, ls="--", lw=1.2, color="0.35")
    ax.set_xlabel("time in the anneal [ps]")
    ax.set_ylim(-1.5, 14.5)
    if k == 0:
        ax.set_ylabel("PE excess over the ruler of the same barostat [meV/atom]")
        ax.legend(fontsize=7.5, loc="center right")
        ax.set_title("(b) 800 K: excess over the ruler", fontsize=10)
    else:
        ax.set_title(f"({'bcd'[k]}) {T} K: excess over the ruler", fontsize=10)
    ax.grid(True, ls="--", color="0.85")
    # caixa dos % liberados (uma linha por família, na cor da família)
    for i, (cor, l) in enumerate(lib_txt):
        ax.text(0.965, 0.052 + i * 0.045, l, transform=ax.transAxes, fontsize=6.5,
                color=cor, ha="right", va="bottom",
                bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none", alpha=0.75))
    ax.text(0.965, 0.190, "released to 1 ns:", transform=ax.transAxes, fontsize=6.8,
            color="0.25", ha="right", va="bottom")

fig.tight_layout()
fig.savefig(SAIDA, dpi=200, bbox_inches="tight")
print("OK:", SAIDA)