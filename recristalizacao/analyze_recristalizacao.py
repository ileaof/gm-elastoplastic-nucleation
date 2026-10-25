# -*- coding: utf-8 -*-
"""Análise das corridas E6 (recristalização por MD com a energia plástica armazenada).

Figuras do PLANO.md, no estilo dos analyze_*.py do Nucleation_MD (matplotlib, Agg):
  fig_E6_estocada.png      E_st por átomo e σ(ε) contra a deformação (Jelinek x Mendelev)
  fig_E6_recozimento.png   PE/átomo e fração FCC (CNA) no tempo, de cada patamar de T
  fig_E6_retratos.png      retratos do fim de cada patamar de T (CNA por cor)
  fig_E6_perfil_delta.png  perfil z da energia por átomo (a espessura δ da zona rica em energia)
Uso: python analyze_recristalizacao.py   (arquivos faltantes são ignorados: sai o que existir)
"""
import glob
import os
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
POTS = ('jelinek', 'mendelev')
COLOR = {'jelinek': '#2a78d6', 'mendelev': '#eb6834'}
TS = 0.002  # ps por passo (units metal)


def read_stored(pot):
    """stored_<pot>.dat: linhas 'STORED eps e_at temp pxy lx' -> (eps, e_at, pxy)."""
    f = os.path.join(HERE, f'stored_{pot}.dat')
    if not os.path.exists(f):
        return None
    eps, eat, pxy = [], [], []
    for line in open(f, encoding='utf-8', errors='replace'):
        p = line.split()
        if p and p[0] == 'STORED':
            eps.append(float(p[1]))
            eat.append(float(p[2]))
            pxy.append(float(p[4]))
    return (np.array(eps), np.array(eat), np.array(pxy)) if eps else None


def thermo_from_log(log):
    """Blocos de termo do log do LAMMPS: (t_ps, temp, pe, pxy) do bloco mais longo."""
    if not os.path.exists(log):
        return None
    raw = open(log, encoding='utf-8', errors='replace').read()
    best = None
    for m in re.finditer(r'^\s*Step\s+Temp\s+Press\s+PotEng[^\n]*\n(.*?)(?:\n\s*\n|\n\s*Loop time)',
                         raw, re.M | re.S):
        rows = [r for r in m.group(1).splitlines() if re.match(r'^\s*\d+\s', r)]
        if len(rows) > (0 if best is None else len(best)):
            best = np.array([[float(x) for x in r.split()[:5]] for r in rows])
    if best is None or best.shape[0] < 3:
        return None
    return best[:, 0] * TS, best[:, 1], best[:, 3], best[:, 4]


def dump_frames(path):
    """Todos os quadros de um dump em formato LAMMPS: lista de (step, frac_FCC, e_med)."""
    if not os.path.exists(path):
        return []
    lines = [l.rstrip('\n') for l in open(path, encoding='utf-8', errors='replace')]
    out, i, cols = [], 0, None
    while i < len(lines):
        if lines[i].startswith('ITEM: TIMESTEP'):
            step, i = int(lines[i + 1]), i + 2
        elif lines[i].startswith('ITEM: NUMBER OF ATOMS'):
            n, i = int(lines[i + 1]), i + 2
        elif lines[i].startswith('ITEM: BOX BOUNDS'):
            i += 4
        elif lines[i].startswith('ITEM: ATOMS'):
            cols = {c: k for k, c in enumerate(lines[i].split()[2:])}
            i += 1
            try:
                arr = np.array([[float(x) for x in lines[i + k].split()] for k in range(n)])
            except (ValueError, IndexError):
                i = len(lines)
                continue
            i += n
            if arr.shape[1] > max(cols.values()):
                out.append((step, float(np.mean(arr[:, cols['c_cnatmp']] == 0)),
                            float(np.mean(arr[:, cols['c_esp']]))))
        else:
            i += 1
    return out


def dump_last_frame(path):
    """Último quadro de um dump: (dict colunas->índice, array id/type/xyz/cna/esp)."""
    fr = dump_frames(path)
    if not fr:
        return None, None
    lines = [l.rstrip('\n') for l in open(path, encoding='utf-8', errors='replace')]
    i, cols, cur = 0, None, None
    step = None
    while i < len(lines):
        if lines[i].startswith('ITEM: TIMESTEP'):
            step = int(lines[i + 1]); i += 2
        elif lines[i].startswith('ITEM: NUMBER OF ATOMS'):
            n = int(lines[i + 1]); i += 2
        elif lines[i].startswith('ITEM: BOX BOUNDS'):
            i += 4
        elif lines[i].startswith('ITEM: ATOMS'):
            cols = {c: k for k, c in enumerate(lines[i].split()[2:])}
            i += 1
            try:
                cur = np.array([[float(x) for x in lines[i + k].split()] for k in range(n)])
            except (ValueError, IndexError):
                cur = None
            i += n
        else:
            i += 1
    return (cols, cur) if cur is not None else (None, None)


def fig_estocada():
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 4.2))
    made = False
    for pot in POTS:
        r = read_stored(pot)
        if r is None:
            continue
        made = True
        eps, eat, pxy = r
        a.plot(eps, eat - eat[0], 'o-', color=COLOR[pot], label=pot)
        b.plot(eps, pxy * 0.1, 'o-', color=COLOR[pot], label=pot)  # bar -> MPa
    if not made:
        return None
    a.axvline(0.14, color='k', ls='-.', lw=1)
    a.set_xlabel('deformação de cisalhamento ε')
    a.set_ylabel('E armazenada Δe/átomo [eV]')
    b.set_xlabel('deformação de cisalhamento ε')
    b.set_ylabel('σ_xy [MPa]')
    a.set_title('Energia plástica armazenada (E6.1)', fontsize=11)
    b.set_title('Tensão × deformação: o escoamento no MD', fontsize=11)
    a.legend()
    for ax in (a, b):
        ax.grid(True, ls='--', color='0.85')
    fig.tight_layout()
    return fig


def fig_recozimento():
    fig, axs = plt.subplots(2, len(POTS), figsize=(6.4 * len(POTS), 6.8), squeeze=False)
    cmap = plt.get_cmap('viridis')
    made = False
    for i, pot in enumerate(POTS):
        logs = sorted(glob.glob(os.path.join(HERE, f'recozer_{pot}_T*.log')))
        for k, log in enumerate(logs):
            T = float(re.search(r'_T([0-9.]+)\.log', log).group(1))
            th = thermo_from_log(log)
            if th is None:
                continue
            made = True
            t_ps, temp, pe, pxy = th
            c = cmap(k / max(len(logs) - 1, 1))
            axs[0][i].plot(t_ps, pe, color=c, lw=1.2, label=f'T = {T:g} K')
            dump = os.path.join(HERE, os.path.basename(log).replace('.log', '.dump')
                                .replace('recozer', 'recozido'))
            fr = dump_frames(dump)
            if fr:
                s, f, e = zip(*fr)
                axs[1][i].plot(np.array(s) * TS, f, '.-', color=c, ms=3, lw=0.8)
        axs[0][i].set_title(f'{pot}: subida 300→T e patamar', fontsize=10)
        axs[0][i].set_xlabel('tempo [ps]')
        axs[0][i].set_ylabel('PE/átomo [eV]')
        axs[1][i].set_xlabel('tempo [ps]')
        axs[1][i].set_ylabel('fração FCC (CNA)')
        axs[1][i].set_ylim(0, 1.05)
        if axs[0][i].get_legend_handles_labels()[0]:
            axs[0][i].legend(fontsize=7)
        for row in axs[:, i]:
            row.grid(True, ls='--', color='0.85')
    fig.tight_layout()
    return fig if made else None


def fig_retratos():
    Ts = [500, 600, 700, 800, 850, 900]
    fig, axs = plt.subplots(len(POTS), len(Ts), figsize=(2.2 * len(Ts), 2.4 * len(POTS)),
                            squeeze=False)
    made = False
    for i, pot in enumerate(POTS):
        for j, T in enumerate(Ts):
            ax = axs[i][j]
            cols, fr = dump_last_frame(os.path.join(HERE, f'recozido_{pot}_T{T}.dump'))
            if fr is None:
                ax.axis('off')
                continue
            made = True
            ax.scatter(fr[:, cols['x']], fr[:, cols['z']], s=1.6,
                       c=fr[:, cols['c_cnatmp']], cmap='coolwarm', vmin=0, vmax=4)
            ax.set_title(f'{pot} {T:g} K', fontsize=8)
            ax.set_xticks([]), ax.set_yticks([])
    if not made:
        return None
    fig.suptitle('CNA do fim de cada patamar: FCC (azul) x defeitos (vermelho)', fontsize=11)
    fig.tight_layout()
    return fig


def fig_perfil_delta():
    fig, ax = plt.subplots(figsize=(8, 4))
    made = False
    for pot in POTS:
        f = os.path.join(HERE, f'perfil_Es_{pot}.dat')
        if not os.path.exists(f):
            continue
        rows = [line.split() for line in open(f, encoding='utf-8', errors='replace')
                if len(line.split()) == 4 and not line.lstrip().startswith('#')]
        if not rows:
            continue
        arr = np.array(rows, float)
        steps = np.unique(arr[:, 0])
        last = arr[arr[:, 0] == steps[-1]]
        ax.plot(last[:, 1], last[:, 3], color=COLOR[pot],
                label=f'{pot} (passo {int(steps[-1])})')
        made = True
    if not made:
        return None
    ax.set_xlabel('z normalizado')
    ax.set_ylabel('E local por átomo [eV]')
    ax.set_title('Perfil de energia: a espessura δ da zona rica em energia', fontsize=11)
    ax.legend()
    ax.grid(True, ls='--', color='0.85')
    fig.tight_layout()
    return fig


def main():
    made = []
    for nome, fn in (('fig_E6_estocada.png', fig_estocada),
                     ('fig_E6_recozimento.png', fig_recozimento),
                     ('fig_E6_retratos.png', fig_retratos),
                     ('fig_E6_perfil_delta.png', fig_perfil_delta)):
        fig = fn()
        if fig is not None:
            fig.savefig(os.path.join(HERE, nome), dpi=150)
            made.append(nome)
    print('figuras:', ', '.join(made) if made else 'nenhuma (sem dados)')


if __name__ == '__main__':
    main()