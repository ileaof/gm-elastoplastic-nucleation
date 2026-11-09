# -*- coding: utf-8 -*-
"""Mapas de grãos/fases por PTM (OVITO) sobre os dumps do E6: cores por fase (FCC/HCP/BCC/outro)
e cores por GRAO (orientacao cristalina) — a recristalizada aparece como grao novo (cor nova).
Uso:  python analyze_graos.py dump1.dump dump2.dump ...   (padrao: os recozidos do acervo)
Saidas: fig_graos_<nome>.png por dump (matplotlib, Agg)."""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ovito.modifiers import PolyhedralTemplateMatchingModifier, StructureIdentificationModifier as _SID
from ovito.io import import_file

CMAP = plt.get_cmap('tab20')
PHASE = {0: ('FCC', '#2a78d6'), 1: ('HCP', '#eb6834'), 2: ('BCC', '#2f8a3a'), 3: ('ICO', '#8a3ffc'), 4: ('other', '0.6')}


def mapa(dump, out_png=None):
    p = import_file(dump)
    data = p.compute(0)
    ptm = PolyhedralTemplateMatchingModifier()
    ptm.structure_types.by_type_name = {0: 'FCC', 1: 'HCP', 2: 'BCC', 3: 'ICO', 4: 'Other'}
    # rodar PTM e granular por orientacao: cluster dos vetores de vizinhanca
    data2 = data.apply(PolyhedralTemplateMatchingModifier())
    pos = data2.particles.positions
    n = len(pos)
    st = data2.particles['Structure Type'][:]
    # cores por grão: cluster de FCC com mesma orientação local (o PTM dá o RMSD e o tipo; a
    # orientação sai da própria geometria: uso o cluster do OVITO sobre o PTM output 'Cluster')
    from ovito.modifiers import ClusterAnalysisModifier
    data3 = data2.apply(ClusterAnalysisModifier(cutoff=4.6, expression="StructureType>=0"))
    cid = data3.particles['Cluster'][:]
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, (label, vals) in zip(axs, (('Fases (PTM)', st), ('Grãos (orientação)', cid))):
        if label.startswith('Fase'):
            cols = np.array([PHASE.get(int(v), ('', '0.4'))[1] for v in vals])
        else:
            uniq = np.unique(vals)
            lut = {u: CMAP((i * 2) % 20) for i, u in enumerate(uniq)}
            cols = np.array([lut[v] for v in vals])
        ax.scatter(pos[:, 0], pos[:, 2], s=2.2, c=cols, lw=0)
        ax.set_title(label)
        ax.set_xticks([])
        ax.set_yticks([])
    if out_png is None:
        out_png = 'fig_graos_' + os.path.basename(dump).replace('.dump', '.png')
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    uniq, cnt = np.unique(st, return_counts=True)
    stat = ', '.join(f'{PHASE[int(u)][0]}={c}' for u, c in zip(uniq, cnt))
    grains = len(np.unique(cid))
    print(f'{os.path.basename(dump)}: {n} átomos | {stat} | {grains} grãos -> {out_png}')
    return out_png


if __name__ == '__main__':
    files = sys.argv[1:] or ['recozido_mendelev_T900.dump', 'recozido_mendelev_T800.dump',
                             'recozido_mendelev_T500.dump', 'deformed_mendelev.dump']
    for f in files:
        f = f if os.path.exists(f) else os.path.join(os.path.dirname(os.path.abspath(__file__)), f)
        if os.path.exists(f):
            try:
                mapa(f)
            except Exception as x:
                print(f, 'ERRO:', x)
        else:
            print('sem arquivo:', f)