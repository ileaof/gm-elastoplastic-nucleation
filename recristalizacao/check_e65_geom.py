# E6.5 — geometria dos vizinhos com mínima imagem exata (fracionária, inclui tilt)
# Perguntas: (a) os átomos "5" têm a casca de 12 vizinhos em ~2,86 Å? (b) onde está a 2ª casca?
# (c) os rótulos CNA formam manchas (estrutura real) ou mistura uniforme (artefato de vibração)?
import os
import numpy as np

E65 = os.path.expanduser("~/e6_rec/reciz_mendelev_T900.dump")

def ler_frame(caminho, buscar_step=None):
    with open(caminho) as f:
        lines = f.read().split("\n")
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: TIMESTEP") and (buscar_step is None or int(lines[i + 1]) == buscar_step):
            natoms = int(lines[i + 3])
            xlo, xhi, xy = [float(v) for v in lines[i + 5].split()]
            ylo, yhi, xz = [float(v) for v in lines[i + 6].split()]
            zlo, zhi, yz = [float(v) for v in lines[i + 7].split()]
            arr = np.array([l.split() for l in lines[i + 9:i + 9 + natoms]], dtype=float)
            cell = np.array([[xhi - xlo, 0.0, 0.0], [xy, yhi - ylo, 0.0], [xz, yz, zhi - zlo]])
            return arr, cell
        i += 1
    raise RuntimeError("frame não achado")

def analisar(arr, cell, tag):
    pos = arr[:, 2:5]
    cna = arr[:, 6].astype(int)
    pe = arr[:, 5]
    inv = np.linalg.inv(cell)
    frac = pos @ inv
    frac -= np.floor(frac)                    # wrap fracionário
    print(f"\n===== {tag} =====")
    print("histograma CNA:", {int(k): int(v) for k, v in zip(*np.unique(cna, return_counts=True))})
    hist_bins = np.linspace(2.0, 4.3, 47)
    print(" rótulo |  n<3.1 |  n<3.6 | d11(Å) | PE méd.")
    for rot in sorted(np.unique(cna)):
        idxs = np.where(cna == rot)[0]
        if len(idxs) == 0:
            continue
        idxs = np.random.default_rng(11).choice(idxs, size=min(300, len(idxs)), replace=False)
        n31, n36, d11, hist = [], [], [], np.zeros(len(hist_bins) - 1)
        for j in idxs:
            dfrac = frac - frac[j]
            dfrac -= np.round(dfrac)              # mínima imagem exata em fracionário
            v = np.sqrt(((dfrac @ cell) ** 2).sum(1))
            v = np.sort(v[v > 1e-2])[:30]
            n31.append(int((v < 3.1).sum())); n36.append(int((v < 3.6).sum()))
            v31 = v[v < 3.1]
            if len(v31) >= 12:
                d11.append(v31[10])
            hist += np.histogram(v, bins=hist_bins)[0]
        n31, n36, d11 = np.array(n31), np.array(n36), np.array(d11)
        print(f"  {rot:4d}   | {n31.mean():6.2f} | {n36.mean():6.2f} | "
              f"{np.median(d11) if len(d11) else float('nan'):6.3f} | {pe[idxs].mean():7.3f}")
        if rot == sorted(np.unique(cna))[0]:
            centros = 0.5 * (hist_bins[:-1] + hist_bins[1:])
            picos = ", ".join(f"{c:.2f}:{v:.0f}" for c, v in zip(centros, hist) if v > hist.max() * 0.08)
            print("  RDF até ~2ª casca (picos):", picos)

    # mistura de rótulos: vizinhos < 3.2 Å — manchas ou mistura uniforme?
    print(" mistura de rótulos nos vizinhos (<3.2 Å):")
    for rot in sorted(np.unique(cna)):
        idxs = np.where(cna == rot)[0]
        if len(idxs) == 0:
            continue
        idxs = np.random.default_rng(11).choice(idxs, size=min(200, len(idxs)), replace=False)
        cont, n_usado = {}, 0
        for j in idxs:
            dfrac = frac - frac[j]
            dfrac -= np.round(dfrac)
            dist = np.sqrt(((dfrac @ cell) ** 2).sum(1))
            viz = np.where((dist > 1e-2) & (dist < 3.2))[0]
            if len(viz) > 16:
                continue                       # ignora casos confusos
            for lb in cna[viz]:
                cont[int(lb)] = cont.get(int(lb), 0) + 1
            n_usado += 1
        tot = sum(cont.values()) or 1
        s = ", ".join(f"{k}: {v / tot * 100:.0f}%" for k, v in sorted(cont.items()))
        print(f"   vizinhos de átomos '{rot}' (n={n_usado}): {s}")

arr, cell = ler_frame(E65, 25000)
analisar(arr, cell, "E6.5 frame 1 (20 ps do recozimento 900 K)")
arr, cell = ler_frame(E65, 500000)
analisar(arr, cell, "E6.5 frame 20 (970 ps do recozimento 900 K)")

for f in ("ctrl_fcc_T300.dump", "ctrl_fcc_T900.dump", "ctrl_fcc_T900_late.dump"):
    p = os.path.expanduser(f"~/e6_rec/{f}")
    if os.path.exists(p) and os.path.getsize(p) > 0:
        try:
            arr, cell = ler_frame(p)   # único frame
            analisar(arr, cell, f"CONTROLE {f}")
        except Exception as e:
            print(f"CONTROLE {f}: erro {e}")