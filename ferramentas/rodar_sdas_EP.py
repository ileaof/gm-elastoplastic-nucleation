"""
Runs the SDAS copies of modelos_elasto-plastos/sdas and saves their figures in sdas/figuras/
(plt.show() of the scripts is replaced by savefig). Usage:  python rodar_sdas_EP.py
"""
import os
import sys
import runpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
SDAS = os.path.join(HERE, '..', 'sdas')
FIG = os.path.join(SDAS, 'figuras')
os.makedirs(FIG, exist_ok=True)
sys.path.insert(0, SDAS)

for f in sorted(os.listdir(SDAS)):
    if not f.endswith('_EP.py'):
        continue
    stem = f[:-3].replace('scriptsdasmodel', '').replace('_zoom_horizontal', '')

    def show(*a, _stem=stem, **k):
        for n in plt.get_fignums():
            plt.figure(n).savefig(os.path.join(FIG, f'{_stem}_fig{n}.png'), dpi=200, bbox_inches='tight')
        plt.close('all')

    plt.show = show
    cwd = os.getcwd()
    os.chdir(SDAS)
    try:
        ns = runpy.run_path(f, run_name='__main__')
        print(f'{f}: Gamma_2nd_EP = {ns["Gibbs_Thomson_het_2nd_EP"]:.4e} m K, '
              f'SDAS_EP = {float(__import__("numpy").exp(ns["z3n_EP"][0])):.2f} t^{ns["z3n_EP"][1]:.3f} um')
    finally:
        os.chdir(cwd)
