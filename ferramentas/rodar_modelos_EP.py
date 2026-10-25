"""
Runs every elasto-plastic copy (sigma_y of its alloy, updated reference, H = 0, loaded surface) in
modelos_elasto-plastos/nucleacao/resultados/<copy>/ and collects the 2nd-order tensor at the mean
experimental gradient (interpolated row of the tables). Writes gamma_2nd_EP_ligas.txt.
Usage:  python rodar_modelos_EP.py
"""
import os
import sys
import csv
import re
import numpy as np
import subprocess
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(__file__))
from aplicar_elastoplastico import MODELS, DST

RES = os.path.join(DST, 'resultados')


def run(item):
    rel, alloy = item
    tag = os.path.basename(rel)[:-3].replace(' ', '_')
    script = os.path.join(DST, tag + '_EP.py')
    short = tag.replace('Nucleation_model_FCC_A1_new_paper_', '').replace('Nucleation_model_', '')[:60]
    wd = os.path.join(RES, short)
    os.makedirs(wd, exist_ok=True)
    env = dict(os.environ, MPLBACKEND='Agg', PYTHONPATH=DST)
    rc = 0
    if '--extract' not in sys.argv:
        with open(os.path.join(wd, 'run.log'), 'w', encoding='utf-8', errors='replace') as log:
            rc = subprocess.run([sys.executable, script], cwd=wd, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=1800).returncode
    G_mean = float(re.search(r'^dTdrj\s*=\s*([-+0-9.eE]+)', open(script, encoding='utf-8').read(), flags=re.M).group(1))
    gt2 = dTdr = None
    # oldest to newest: the last table read wins (older runs of other energy modes stay in the folder)
    for f in sorted(os.listdir(wd), key=lambda n: os.path.getmtime(os.path.join(wd, n))):
        if f.endswith('.csv') and (f.startswith('tablefull') or f.startswith('tabela')):
            rows = list(csv.reader(open(os.path.join(wd, f)), delimiter=';'))
            head = [h.strip() for h in rows[0]]
            if 'GT_het_2nd' not in head:
                continue
            interp = [r for r in rows[1:] if r and r[0].strip().endswith('.5')]
            if interp:
                gt2 = float(interp[-1][head.index('GT_het_2nd')])
                dTdr = float(interp[-1][head.index('dTdr')])
            else:
                # no interpolated row in this script: interpolate at its mean gradient dTdrj
                data = [r for r in rows[1:] if r and r[0].strip().isdigit() and int(r[0]) > 0]
                names = [h for h in head if h]
                k2 = head.index('GT_het_2nd')
                if len(data[0]) == len(names) + 1:
                    # header of the original 'tabela' has one label less than the columns ('GT_het' missing):
                    # the column labelled GT_het_2nd holds GT_het and GT_het_2nd is the last, unlabelled one
                    k2 = len(names)
                G = np.array([float(r[head.index('dTdr')]) for r in data])
                g = np.array([float(r[k2]) for r in data])
                gt2, dTdr = float(np.interp(G_mean, G, g)), G_mean
    sweep = [l for l in open(os.path.join(wd, 'run.log'), encoding='utf-8', errors='replace') if l.startswith('Sweep:')]
    return alloy, short, rc, dTdr, gt2, sweep[0].strip() if sweep else ''


lines = ['# copy  alloy  grad_T_mean[K/m]  Gamma_2nd_EP[m K]  (rodar_modelos_EP.py)']
with ThreadPoolExecutor(max_workers=6) as ex:
    for alloy, short, rc, dTdr, gt2, sweep in ex.map(run, MODELS.items()):
        print(f'{short:<62} {alloy:<9} exit {rc}  grad T = {dTdr}  Gamma_2nd_EP = {gt2}  | {sweep}')
        lines.append(f'{short} {alloy} {dTdr} {gt2}')
open(os.path.join(os.path.dirname(__file__), 'gamma_2nd_EP_ligas.txt'), 'w').write('\n'.join(lines) + '\n')
