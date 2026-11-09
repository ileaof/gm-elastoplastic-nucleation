# E6.5 — diagnóstico do CNA: os átomos "5" têm casca de 12 vizinhos FCC?
# Pega o 1º e o 20º frame, amostra átomos por rótulo CNA, mede os vizinhos até 3,6 Å,
# e testa a ordem estrutural (distâncias da casca, ângulo/coordenação).
import os
import numpy as np

DUMP = os.path.expanduser("~/e6_rec/reciz_mendelev_T900.dump")
RC = 3.6   # um pouco acima da 1ª casca (2,86), abaixo da 2ª (4,045)

os.makedirs(os.path.expanduser("~/e6_rec/e65_analise"), exist_ok=True)

def ler_frame(buscar_step):
    with open(DUMP) as f:
        lines = f.read().split("\n")
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: TIMESTEP") and int(lines[i + 1]) == buscar_step:
            natoms = int(lines[i + 3])
            arr = np.array([l.split() for l in lines[i + 9:i + 9 + natoms] if l.strip()], dtype=float)
            return arr
        i += 1
    raise RuntimeError("frame não achado")

def vizinhos(arr, amostra=400, seed=7):
    pos = arr[:, 2:5]
    ids = arr[:, 0].astype(int)
    cna = arr[:, 6].astype(int)
    rng = np.random.default_rng(seed)
    por_rotulo = {}
    for rot in np.unique(cna):
        idx_rot = np.where(cna == rot)[0]
        if len(idx_rot) == 0:
            continue
        escolha = rng.choice(idx_rot, size=min(amostra, len(idx_rot)), replace=False)
        cnts, n11, resto = [], [], 0.0
        for j in escolha:
            d = np.linalg.norm(pos - pos[j], axis=1)
            vasc = (d > 1e-3) & (d < RC)
            n = vasc.sum()
            dd = np.sort(d[vasc])
            cnts.append(n)
            if len(dd) >= 11:
                n11.append(dd[10])
        por_rotulo[int(rot)] = (len(escolha), np.mean(cnts), np.std(cnts),
                                np.min(cnts), np.max(cnts), np.mean(n11) if n11 else np.nan)
    return por_rotulo

for step in (25000, 500000):
    arr = ler_frame(step)
    print(f"\n=== frame tstep {step} ===")
    print("histograma CNA do frame inteiro:",
          {int(k): int(v) for k, v in zip(*np.unique(arr[:, 6].astype(int), return_counts=True))})
    res = vizinhos(arr)
    print("rótulo | n_amostrados | vizinhos méd. | σ | mín | máx | 11º vizinho (Å)")
    for rot, (n, m, s, lo, hi, v11) in sorted(res.items()):
        print(f"{rot:6d} | {n:6d} | {m:6.2f} | {s:4.2f} | {lo:4.0f} | {hi:4.0f} | {v11:6.3f}")
    # padrão espacial: fração de '5' por fatia em x (planos? faixas? aleatório?)
    cna = arr[:, 6].astype(int)
    fatias = np.linspace(arr[:, 2].min(), arr[:, 2].max(), 11)
    prof = []
    for a, b in zip(fatias[:-1], fatias[1:]):
        m = (arr[:, 2] >= a) & (arr[:, 2] < b)
        prof.append((cna[m] == 5).mean() if m.any() else np.nan)
    print("fração '5' por fatia em x:", " ".join(f"{v:.2f}" for v in prof))