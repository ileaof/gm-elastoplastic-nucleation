# Fig 10 do paper, versão EN (eixos em inglês, o mesmo gráfico): onset do decaimento catastrófico
# por temperatura/semente, DO PARSE DA PRÓPRIA TABELA 9 (g3_events): o dpe[0] do início do run
# é saltado (artefato da linha de partida — o mesmo salto que dava "0 ps" no s11-775 e "1 ps/2 ps" no longo).
import glob
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
SAIDA = r"C:\Users\ileao\OneDrive\Documentos\gm-elastoplastic-nucleation\paper\Fig10_E6_onset.png"
EP6 = os.path.join(os.path.dirname(AQUI), "recristalizacao")

def eventos():
    out = []
    for log in sorted(glob.glob(os.path.join(EP6, "gpu2", "recozer_jelinek_T*_s*.log")) +
                      glob.glob(os.path.join(EP6, "recozer_jelinek_T*_s*.log"))):
        mm = re.match(r"recozer_jelinek_T([0-9.]+)_s([0-9]+)\.log", os.path.basename(log))
        if not mm:
            continue
        T, S = float(mm.group(1)), int(mm.group(2))
        if [(r[0], r[1]) for r in out] and (T, S) in [(float(r[0]), int(r[1])) for r in out]:
            continue
        best = None
        for m in re.finditer(r"^\s*Step\s+Temp\s+Press\s+PotEng[^\n]*\n(.*?)(?:\n\s*\n|\n\s*Loop time|\Z)",
                             open(log, encoding="utf-8", errors="replace").read(), re.M | re.S):
            rr = [r.split() for r in m.group(1).splitlines()
                  if re.match(r"^\s*\d+\s", r) and len(r.split()) == 6]
            if rr and (best is None or len(rr) > len(best)):
                best = np.array([[float(x) for x in r] for r in rr])
        if best is None:
            continue
        step, temp, pe = best[:, 0], best[:, 1], best[:, 3]
        dpe = np.diff(pe)[1:]          # salta o dpe da linha de partida (aquecimento, não evento)
        i = int(np.argmax(dpe)) + 1
        raw = open(log, encoding="utf-8", errors="replace").read()
        lost = re.search(r"Lost atoms: original (\d+) current (\d+)", raw)
        n_lost = int(lost.group(1)) - int(lost.group(2)) if lost else 0
        if n_lost:
            out.append([T, S, step[i] * 0.002, temp.max()])
        else:
            out.append([T, S, None, None])
    return out

ev = eventos()
print("eventos parseados:", ev)

fig, ax = plt.subplots(figsize=(9.6, 4.6))
ax.axhline(100, color="k", lw=1.2)
LIM_T = (735, 815)
for T, S, t_ev, tmáx in ev:
    if t_ev is None:
        ax.plot(T, 100, "s", ms=10, color="#3a6ea5", zorder=3)
        ax.annotate(f"s{S} clean", (T, 100), textcoords="offset points", xytext=(0, 12),
                    fontsize=8.5, color="#3a6ea5", ha="center", zorder=5)
    else:
        ax.plot(T, t_ev, "o", ms=9, color="#d95f02", zorder=3)
        dy = 6 if t_ev < 50 else 4
        ax.annotate(f"s{S}", (T, t_ev), textcoords="offset points",
                    xytext=(14, dy), fontsize=8.5, color="#a8501f")
ax.set_xlabel("annealing temperature [K]")
ax.set_ylabel("onset time of the decay [ps]")
ax.set_xlim(*LIM_T)
ax.set_ylim(12, 106)
ax.set_title("The onset of the catastrophic decay at 750-800 K: 6 of 9 anneals decayed "
             "within 17-66 ps", fontsize=11)
ax.grid(True, ls="--", color="0.85")
fig.tight_layout()
fig.savefig(SAIDA, dpi=200)
print("OK:", SAIDA)