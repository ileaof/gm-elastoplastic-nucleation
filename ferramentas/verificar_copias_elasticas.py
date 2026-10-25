"""
Check that every elasto-plastic copy, run with GM_SURFACE_MODEL=elastic, reproduces its original.
Each original and each copy runs in its own temporary folder; the CSV tables are compared row by row
(rows with the same first column, i.e. the same r step, and the interpolated row at the experimental
gradient). Usage:  python verificar_copias_elasticas.py <scratch folder>
"""
import os
import sys
import csv
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(__file__))
from aplicar_elastoplastico import MODELS, SRC, DST, HELPERS

scratch = os.path.abspath(sys.argv[1])
os.makedirs(scratch, exist_ok=True)


def run(script_path, workdir, env_extra):
    os.makedirs(workdir, exist_ok=True)
    for h in HELPERS:
        shutil.copy2(os.path.join(SRC, h), workdir)
    shutil.copy2(os.path.join(DST, 'gm_elastoplastic.py'), workdir)
    shutil.copy2(script_path, os.path.join(workdir, 'model.py'))
    env = dict(os.environ, MPLBACKEND='Agg', **env_extra)
    with open(os.path.join(workdir, 'run.log'), 'w', encoding='utf-8', errors='replace') as log:
        p = subprocess.run([sys.executable, 'model.py'], cwd=workdir, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=1200)
    return p.returncode


def tables(workdir):
    out = {}
    for f in os.listdir(workdir):
        if f.endswith('.csv') and (f.startswith('tablefull') or f.startswith('tabela')):
            kind = 'tablefull' if f.startswith('tablefull') else 'tabela'
            rows = list(csv.reader(open(os.path.join(workdir, f)), delimiter=';'))
            out[kind] = {r[0].strip(): r for r in rows[1:] if r and r[0].strip()}
    return out


def compare(a, b):
    worst, n = 0.0, 0
    for key, ra in a.items():
        rb = b.get(key)
        if rb is None:
            continue
        n += 1
        for x, y in zip(ra, rb):
            try:
                x, y = float(x), float(y)
            except ValueError:
                continue
            if x != x and y != y:
                continue
            worst = max(worst, abs(x - y) / max(abs(x), 1e-30))
    return n, worst


def job(item):
    idx, (rel, alloy) = item
    tag = os.path.basename(rel)[:-3].replace(' ', '_')
    copy = os.path.join(DST, tag + '_EP.py')
    # short folder names: Windows paths are limited to 260 characters
    w_o, w_c = os.path.join(scratch, f'm{idx:02d}', 'o'), os.path.join(scratch, f'm{idx:02d}', 'c')
    rc_o = run(os.path.join(SRC, rel), w_o, {})
    rc_c = run(copy, w_c, {'GM_SURFACE_MODEL': 'elastic'})
    to, tc = tables(w_o), tables(w_c)
    res = []
    for kind in to:
        if kind in tc:
            n, worst = compare(to[kind], tc[kind])
            res.append(f'{kind}: {n} rows, max rel diff {worst:.2e}')
    return tag, rc_o, rc_c, res


with ThreadPoolExecutor(max_workers=6) as ex:
    for tag, rc_o, rc_c, res in ex.map(job, enumerate(MODELS.items())):
        print(f'{tag[:80]:<80} exit orig {rc_o} copy {rc_c} | ' + ' | '.join(res))
