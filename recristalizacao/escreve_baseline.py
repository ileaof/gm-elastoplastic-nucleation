# Grava o baseline dos controles FCC puros (62 500 átomos) em
# <repo>/recristalizacao/e65_serie26/baseline_puro.dat para o build do paper.
import os
import numpy as np

DEST = r"C:\Users\ileao\OneDrive\Documentos\gm-elastoplastic-nucleation\recristalizacao\e65_serie26"
AQUI = os.path.dirname(os.path.abspath(__file__))
CONTROLES = {800: ("ctrl_fcc_T800_hold0.dump", "ctrl_fcc_T800_late.dump"),
             850: ("ctrl_fcc_T850_hold0.dump", "ctrl_fcc_T850_late.dump"),
             900: ("ctrl_fcc_T900.dump", "ctrl_fcc_T900_late.dump")}
valores = {}
for T, arqs in CONTROLES.items():
    f5, pe = [], []
    for a in arqs:
        lines = open(os.path.join(AQUI, a)).read().split("\n")
        n = int(lines[3])
        arr = np.array([l.split() for l in lines[9:9 + n]], dtype=float)
        m = arr[:, 6].astype(int) == 5
        f5.append(m.mean())
        pe.append(arr[:, 5].mean())
    valores[T] = (float(np.mean(f5)), float(np.mean(pe)))
    print(T, f"{np.mean(f5) * 100:.2f}%", f"{np.mean(pe):.4f} eV/atom")

os.makedirs(DEST, exist_ok=True)
with open(os.path.join(DEST, "baseline_puro.dat"), "w") as g:
    g.write("# FCC puro, 62 500 atoms, same computes: CNA label-5 fraction and PE per atom at the\n"
            "# anneal temperature (the per-temperature baseline of the hot-CNA artifact).\n")
    g.write("# T[K] f5_pure PE_pure[eV/atom]\n")
    for T, (f5, pe) in sorted(valores.items()):
        g.write(f"{T} {f5:.6f} {pe:.6f}\n")
print("gravado em", os.path.join(DEST, "baseline_puro.dat"))