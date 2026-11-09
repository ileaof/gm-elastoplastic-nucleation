# E6.5 série — resumo dos controles (f5 e PE baselines por T, hold0 e late reais)
# e do Pxy/PE dos holds 850/800 do rec_gpu_T{850,800}.log
import numpy as np

N = 62500
VOL = 1.034e-24  # m3

def dump(path):
    lines = open(path).read().split("\n")
    n = int(lines[3])
    arr = np.array([l.split() for l in lines[9:9 + n]], dtype=float)
    return arr

print("=== controles FCC puro: baseline por temperatura ===")
print(f"{'T':>6} {'f5(%)':>7} {'PE_media':>9} {'PE_fcc':>9} {'PE_5':>9}")
bases = {}
for T in (800, 850, 900):
    if T == 900:
        paths = {"hold0": "ctrl_fcc_T900.dump", "late": "ctrl_fcc_T900_late.dump"}
    else:
        paths = {"hold0": f"ctrl_fcc_T{T}_hold0.dump", "late": f"ctrl_fcc_T{T}_late.dump"}
    for tag, p in paths.items():
        try:
            a = dump(p)
        except Exception as e:
            print(f"  {T} {tag}: ERRO {e}")
            continue
        cna = a[:, 6].astype(int)
        pe = a[:, 5]
        m5 = cna == 5
        print(f"{T:4d}K  {tag:5} {m5.mean() * 100:6.2f}  {pe.mean():9.4f} "
              f"{pe[~m5].mean():9.4f} {pe[m5].mean():9.4f}")
        bases[(T, tag)] = (m5.mean(), pe.mean())
for T in (800, 850, 900):
    f5 = 100 * (bases.get((T, "hold0"), (np.nan,))[0] + bases.get((T, "late"), (np.nan,))[0]) / 2
    pe = (bases.get((T, "hold0"), (0, np.nan))[1] + bases.get((T, "late"), (0, np.nan))[1]) / 2
    print(f"baseline T={T}: f5 = {f5:.2f}%, PE = {pe:.4f} eV/átomo")

print("\n=== E6.5 holds: Pxy/PE do thermo (rec_gpu_T*.log) ===")
for T in (900, 850, 800):
    nome = "rec_gpu.log" if T == 900 else f"rec_gpu_T{T}.log"
    try:
        rows = []
        for ln in open(nome, errors="ignore"):
            p = ln.split()
            if len(p) == 6:
                try:
                    rows.append([float(x) for x in p])
                except ValueError:
                    pass
        rows = np.array(rows)
        hold = rows[rows[:, 0] >= 15000]
        print(f"T={T}: Pxy médio {hold[:, 4].mean() / 1000:7.1f} kbar (±{hold[:, 4].std() / 1000:.1f}); "
              f"PE médio {hold[:, 3].mean() / N:9.4f} eV/átomo; {len(hold)} amostras thermo")
    except FileNotFoundError:
        print(f"T={T}: {nome} não achado")

print("\n=== excesso e ΔPE (E6.5 - puro) por T ===")
# E6.5: do CSV dos analisadores (média dos frames)
for T in (800, 850, 900):
    try:
        csvdir = "e65_analise" if T == 900 else f"e65_analise_T{T}"
        c = np.genfromtxt(f"{csvdir}/e65_resumo.csv", delimiter=",", names=True)
        f5_e65 = c["f_diamante"].mean() * 100
        pe_e65 = c["PE_medio"].mean()
        b = (bases.get((T, "hold0")), bases.get((T, "late")))
        f5_puro = np.nanmean([x[0] for x in b]) * 100
        pe_puro = np.nanmean([x[1] for x in b])
        d = pe_e65 - pe_puro
        dens = d * N * 1.602176634e-19 / VOL
        print(f"T={T}: f5 bruto {f5_e65:5.1f}% | excesso {f5_e65 - f5_puro:5.1f}% | "
              f"PE {pe_e65:.4f} vs puro {pe_puro:.4f} | Δ {d:+.4f} eV/átomo ≈ {dens:.2e} J/m³")
    except Exception as e:
        print(f"T={T}: erro {e}")