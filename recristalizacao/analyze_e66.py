# E6.6 — análise completa na caixa: descarregamentos (cruzamento + hold NPT), réguas
# NPT/NVT e corridas largas (γ=0.20/0.30 NVT). Saídas CSV em ~/e6_rec/e66_analise/.
import glob
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.expanduser("~/e6_rec")
SAIDA = os.path.join(AQUI, "e66_analise")
N = 62500
AREA = (101.125e-10) ** 2   # plano y-normal [m2] (Lx×Lz, como o E6.4)
J = 1.602177e-19
DT = 0.002

os.makedirs(SAIDA, exist_ok=True)

def thermo_serie(caminho):
    rows = []
    for ln in open(caminho, errors="replace"):
        p = ln.split()
        if len(p) == 6:
            try:
                rows.append([float(x) for x in p])
            except ValueError:
                pass
    return np.array(rows)

def dump_frame(caminho, buscar_step=None):
    with open(caminho) as f:
        lines = f.read().split("\n")
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: TIMESTEP") and (buscar_step is None or int(lines[i + 1]) == buscar_step):
            n = int(lines[i + 3])
            arr = np.array([l.split() for l in lines[i + 9:i + 9 + n]], dtype=float)
            if arr.shape[0] == 0:
                return None
            return arr
        i += 1
    return None

def dump_serie(caminho, base=None):
    """Cada frame: (t_ps desde base, f5, PE média, PE fcc, PE nao-fcc)."""
    with open(caminho) as f:
        lines = f.read().split("\n")
    res, i = [], 0
    while i < len(lines):
        if lines[i].startswith("ITEM: TIMESTEP"):
            step = int(lines[i + 1])
            n = int(lines[i + 3])
            arr = np.array([l.split() for l in lines[i + 9:i + 9 + n]], dtype=float)
            cna = arr[:, 6].astype(int)
            m5 = cna == 5
            mf = cna == 1
            ref = 15000 if base is None else base
            res.append([max(0.0, (step - ref) * DT), m5.mean(), arr[:, 5].mean(),
                        arr[mf, 5].mean() if mf.any() else np.nan,
                        arr[~m5, 5].mean() if m5.any() else np.nan])
            i += 9 + n
        else:
            i += 1
    return np.array(res)

def lost_atoms(caminho):
    m = re.search(r"Lost atoms: original (\d+) current (\d+)", open(caminho, errors="replace").read())
    return int(m.group(1)) - int(m.group(2)) if m else 0

def hold_da(caminho, passo_ini):
    s = thermo_serie(caminho)
    h = s[s[:, 0] >= passo_ini]
    return s, h

print("========== 1) RÉGUAS ==========")
regua = {}
# NVT: controles puros já medidos (800/850/900 ontem; 500/600/700 hoje)
for T, arqs in ((500, ("ctrl_fcc_T500_hold0.dump", "ctrl_fcc_T500_late.dump")),
                (600, ("ctrl_fcc_T600_hold0.dump", "ctrl_fcc_T600_late.dump")),
                (700, ("ctrl_fcc_T700_hold0.dump", "ctrl_fcc_T700_late.dump")),
                (800, ("ctrl_fcc_T800_hold0.dump", "ctrl_fcc_T800_late.dump")),
                (850, ("ctrl_fcc_T850_hold0.dump", "ctrl_fcc_T850_late.dump")),
                (900, ("ctrl_fcc_T900.dump", "ctrl_fcc_T900_late.dump"))):
    f5, pe = [], []
    for a in arqs:
        arr = dump_frame(os.path.join(AQUI, a))
        if arr is None:
            continue
        m5 = arr[:, 6].astype(int) == 5
        f5.append(m5.mean()); pe.append(arr[:, 5].mean())
    if pe:
        regua[("NVT", T)] = (float(np.mean(f5)), float(np.mean(pe)))
# NPT: controles puros de hoje
for T in (800, 850, 900):
    _, h = hold_da(os.path.join(AQUI, f"ctrl_npt_T{T}.log"), 0)
    fim = h[h[:, 0] >= 20000]     # hold após rampa de 5000 (5000+1+20000+1)
    f5, pe = [], []
    for tag in ("hold0", "late"):
        arr = dump_frame(os.path.join(AQUI, f"ctrl_npt_T{T}_{tag}.dump"))
        if arr is not None:
            m5 = arr[:, 6].astype(int) == 5
            f5.append(m5.mean()); pe.append(arr[:, 5].mean())
    regua[("NPT", T)] = (float(np.mean(f5)) if f5 else np.nan, float(np.mean(pe)) if pe else np.nan,
                         float(fim[:, 2].mean()) / 1000 if len(fim) else np.nan,
                         float(fim[:, 3].mean() / N) if len(fim) else np.nan)

with open(os.path.join(SAIDA, "reguas.csv"), "w") as g:
    g.write("tipo,T,f5_pure,PE_pure,press_pure[kbar],PE_thermo\n")
    for (tipo, T), v in regua.items():
        if tipo == "NPT":
            g.write(f"{tipo},{T},{v[0]:.6f},{v[1]:.6f},{v[2]:.2f},{v[3]:.6f}\n")
            print(f"régua NPT {T} K: f5 {v[0]*100:5.2f}%  PE(dump) {v[1]:.4f}  "
                  f"PE(thermo) {v[3]:.4f}  press {v[2]:+.1f} kbar")
        else:
            g.write(f"{tipo},{T},{v[0]:.6f},{v[1]:.6f},,\n")
            print(f"régua NVT {T} K: f5 {v[0]*100:5.2f}%  PE {v[1]:.4f}")

print("\n========== 2) DESCARREGAMENTO + NPT ==========")
with open(os.path.join(SAIDA, "descarrega.csv"), "w") as g:
    g.write("tag,T,seed,gam,eps_cross,PE_cross,resto_eV_at,resto_Jm2,"
            "PE_hold_ini,PE_hold_fim,dPE_fim,p_evento,Pxy_med[kbar],press_med[kbar],vol_ini,vol_fim,perdeu\n")
    for logp in sorted(glob.glob(os.path.join(AQUI, "desca_npt*.log"))):
        tag = os.path.basename(logp).replace(".log", "")
        mm = re.match(r"desca_npt(30)?_T([0-9]+)(?:_s([0-9]+))?", tag)
        if not mm:
            continue
        gam = "0.30" if mm.group(1) else "0.10"
        T = int(mm.group(2)); seed = mm.group(3) or "101"
        bruto = open(logp, errors="replace").read()
        cr = re.search(r"DESCA_CROSS step=(\d+) pxy=([-\d.]+) PE=(-[\d.]+)", bruto)
        passo_cross = int(cr.group(1)); pe_cross = float(cr.group(3)) / N
        fim_cis = 5000 + int(round(50000 * float(gam)))
        eps_cross = float(gam) - 2e-5 * (passo_cross - fim_cis)
        rows = thermo_serie(logp)
        pe_eq = rows[(rows[:, 0] >= 4999) & (rows[:, 0] <= 5010)][:, 3].mean() / N   # eq 300 K pre-cis
        resto_at = pe_cross - pe_eq
        resto_Jm2 = resto_at * N * J / AREA
        ini_hold = passo_cross + 7000
        _, h = hold_da(logp, ini_hold)
        if not len(h):
            continue
        pe_ini, pe_fim = h[0, 3] / N, h[-1, 3] / N
        d_evento = 0.0
        if len(h) > 1:
            d_evento = float(np.max(np.abs(np.diff(h[:, 3]))) / N)
        perdeu = lost_atoms(logp)
        g.write(f"{tag},{T},{seed},{gam},{eps_cross:.4f},{pe_cross:.5f},{resto_at:.5f},{resto_Jm2:.3f},"
                f"{pe_ini:.5f},{pe_fim:.5f},{pe_fim - pe_ini:+.5f},{d_evento:.4f},"
                f"{h[:, 4].mean() / 1000:.1f},{h[:, 2].mean() / 1000:.2f},"
                f"{h[0, 5]:.1f},{h[-1, 5]:.1f},{perdeu}\n")
        print(f"{tag}: eps_cruz {eps_cross:.3f} (γ={gam}); resto {resto_at:+.5f} eV/at = "
              f"{resto_Jm2:.2f} J/m² | hold PE {pe_ini:.4f}→{pe_fim:.4f} | maior Δpasso {d_evento:.3f} | "
              f"Pxy {h[:, 4].mean() / 1000:.0f} kbar | vol {h[0, 5]:.0f}→{h[-1, 5]:.0f} | perdeu {perdeu}")
        # serie de frames do hold (dumps)
        dmp = os.path.join(AQUI, f"reciz_descarnpt_T{T}_g{gam}_s{seed}.dump")
        if os.path.exists(dmp) and os.path.getsize(dmp):
            fr = dump_serie(dmp, base=passo_cross + 7000)
            np.savetxt(os.path.join(SAIDA, f"frames_descarnpt_T{T}_g{gam}_s{seed}.csv"), fr,
                       delimiter=",", header="t_ps,f5,PE,PE_fcc,PE_nao_fcc", comments="")

print("\n========== 3) LARGAS (γ carregado, NVT) ==========")
with open(os.path.join(SAIDA, "largo.csv"), "w") as g:
    g.write("tag,T,gam,PE_hold_ini,PE_hold_fim,dPE_fim,p_evento,Pxy_med[kbar],f5_ini,f5_fim,perdeu\n")
    for T in (800, 850, 900):
        for gam in ("0.20", "0.30"):
            logp = os.path.join(AQUI, f"rec_largo{int(round(float(gam) * 100))}_T{T}.log")
            if not os.path.exists(logp):
                print(f"largo T={T} γ={gam}: log não achado ({logp})")
                continue
            _, h = hold_da(logp, 5000 + int(round(50000 * float(gam))) + 5000)
            if not len(h):
                continue
            d_evento = float(np.max(np.abs(np.diff(h[:, 3]))) / N) if len(h) > 1 else 0.0
            dmp = os.path.join(AQUI, f"reciz_largo_T{T}_g{gam}.dump")
            fr = dump_serie(dmp, base=10000 + int(round(50000 * float(gam)))) \
                if os.path.exists(dmp) and os.path.getsize(dmp) else None
            if fr is not None:
                np.savetxt(os.path.join(SAIDA, f"frames_largo_T{T}_g{gam}.csv"), fr, delimiter=",",
                           header="t_ps,f5,PE,PE_fcc,PE_nao_fcc", comments="")
            f5_ini = fr[0, 1] if fr is not None else np.nan
            f5_fim = fr[-1, 1] if fr is not None else np.nan
            perdeu = lost_atoms(logp)
            g.write(f"largo{T}_g{gam},{T},{gam},{h[0, 3] / N:.5f},{h[-1, 3] / N:.5f},"
                    f"{h[-1, 3] / N - h[0, 3] / N:+.5f},{d_evento:.4f},{h[:, 4].mean() / 1000:.1f},"
                    f"{f5_ini:.4f},{f5_fim:.4f},{perdeu}\n")
            print(f"largo T={T} γ={gam}: PE {h[0, 3] / N:.4f}→{h[-1, 3] / N:.4f} | "
                  f"f5 {f5_ini:.3f}→{f5_fim:.3f} | Pxy {h[:, 4].mean() / 1000:.0f} kbar | perdeu {perdeu}")

# régua NVT e os E6.5 (γ=0.10 carregado NVT) já calculados antes: e65_resumo_T800/T850/T900
# estão no repo; para a comparação direta, extraio aqui as séries de frames do E6.5:
for T in (800, 850, 900):
    dmp = os.path.join(AQUI, f"reciz_mendelev_T{T}.dump")
    if os.path.exists(dmp) and os.path.getsize(dmp):
        fr = dump_serie(dmp)
        np.savetxt(os.path.join(SAIDA, f"frames_e65_T{T}_g0.10.csv"), fr, delimiter=",",
                   header="t_ps,f5,PE,PE_fcc,PE_nao_fcc", comments="")
        print(f"E6.5 T={T} γ=0.10 frames: {len(fr)}")

print("\nOK — saídas em", SAIDA)