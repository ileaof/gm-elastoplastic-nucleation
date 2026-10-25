"""
Texts of the nucleation figures inside the axes (elasto-plastic copies, modelos_elasto-plastos).

The texts of the original scripts (alloy name, N/V, E_Fermi, ...) are placed in data coordinates chosen
for the original sweep; with the elasto-plastic sweep some of them fall outside the axes. texts_inside_all()
moves every text that leaves its axes or covers the legend into the corner of the axes that overlaps least the legend and
the curves, keeping their order; a legend outside the axes goes back inside (loc='best').
Other texts are not moved.
"""
import numpy as np
import matplotlib.pyplot as plt


def _box(artist, renderer, inv):
    b = artist.get_window_extent(renderer).transformed(inv)
    return (b.x0, b.y0, b.x1, b.y1)


def _leaves(b, tol=0.005):
    """Box (axes coordinates) not contained in the axes."""
    return b[0] < -tol or b[1] < -tol or b[2] > 1 + tol or b[3] > 1 + tol


def _misplaced(ax, tx, renderer, inv, lb):
    """Text in data coordinates that leaves the axes or covers the legend."""
    if tx.get_transform() != ax.transData:
        return False
    b = _box(tx, renderer, inv)
    if _leaves(b):
        return True
    return lb is not None and _overlap(b, lb) > 0.0


def _overlap(a, b):
    w = min(a[2], b[2]) - max(a[0], b[0])
    h = min(a[3], b[3]) - max(a[1], b[1])
    return max(w, 0.0) * max(h, 0.0)


def texts_inside(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in fig.axes:
        inv = ax.transAxes.inverted()
        leg = ax.get_legend()
        if leg is not None and _leaves(_box(leg, renderer, inv)):   # legend outside the axes
            leg.set_bbox_to_anchor(None)
            leg.set_loc('best')
            fig.canvas.draw()
            size = max(t.get_fontsize() for t in leg.get_texts())
            while _leaves(_box(leg, renderer, inv)) and size > 7:   # larger than the axes: smaller font
                size -= 1
                for t in leg.get_texts():
                    t.set_fontsize(min(t.get_fontsize(), size))
                fig.canvas.draw()
        lb = _box(leg, renderer, inv) if leg is not None else None
        outs = [tx for tx in ax.texts if _misplaced(ax, tx, renderer, inv, lb)]
        if not outs:
            continue
        outs.sort(key=lambda tx: -tx.get_position()[1] * (1 if ax.get_ylim()[1] > ax.get_ylim()[0] else -1))
        sizes = []
        for tx in outs:
            tx.set_transform(ax.transAxes)
            tx.set_position((0.0, 0.0))
            tx.set_ha('left')
            tx.set_va('top')
            tx.set_fontsize(min(tx.get_fontsize(), 12))
            bb = tx.get_window_extent(renderer).transformed(inv)
            sizes.append((bb.width, bb.height))
        gap = 0.012
        W = max(w for w, _ in sizes)
        H = sum(h for _, h in sizes) + gap * (len(sizes) - 1)
        pts = []
        for ln in ax.get_lines():
            xy = np.column_stack(ln.get_data()).astype(float)
            if len(xy):
                pts.append(inv.transform(ax.transData.transform(xy)))
        pts = np.vstack(pts) if pts else np.zeros((0, 2))
        best = None
        for x0 in (0.03, max(0.03, 0.97 - W)):
            for ytop in (0.97, 0.03 + H):
                rect = (x0, ytop - H, x0 + W, ytop)
                score = 1000 * _overlap(rect, lb) if lb else 0.0
                if len(pts):
                    score += np.mean((pts[:, 0] >= rect[0]) & (pts[:, 0] <= rect[2]) &
                                     (pts[:, 1] >= rect[1]) & (pts[:, 1] <= rect[3]))
                if best is None or score < best[0]:
                    best = (score, x0, ytop)
        _, x0, y = best
        for tx, (_, h) in zip(outs, sizes):
            tx.set_position((x0, y))
            y -= h + gap


def elastic_rate(script):
    """Nucleation rate of the same script with the elastic surface (run in a separate process, GM_DUMP_RATE),
    for the comparison in the rate figure. Returns a dict of arrays (DT, G, Ihom, Ihet) or None."""
    import os
    import sys
    import subprocess
    import tempfile
    with tempfile.TemporaryDirectory() as wd:
        out = os.path.join(wd, 'rate_elastic.npz')
        env = dict(os.environ, GM_SURFACE_MODEL='elastic', GM_DUMP_RATE=out, MPLBACKEND='Agg')
        try:
            subprocess.run([sys.executable, os.path.abspath(script)], cwd=wd, env=env, timeout=1800,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            with np.load(out) as d:   # closed before the temporary folder is removed (Windows)
                return {k: np.array(d[k]) for k in d.files}
        except Exception as e:
            print(f'figure_layout: elastic rate not available ({type(e).__name__}: {e})')
            return None


def vs_gradient_figure(G, y, ylabel, title, ep, y_el=None, G_el=None, G_yield=None, G_mean=None, G_mark=1e5,
                       log_y=True, unit='', extra=(), main_label=None):
    """One quantity of the nucleation sweep against the thermal gradient (whole range of the sweep): elasto-plastic
    (or elastic) curve of the script, elastic curve of the same script, yield (elastic range shaded), mean experimental
    gradient and the values at G_mark."""
    sfx = ', elasto-plastic surface' if ep else ', elastic surface'
    fig, ax = plt.subplots(1)
    ax.set_title(title, {'color': 'black', 'fontsize': 13})
    if ep and G_yield:
        ax.axvspan(min(G[0], G_el[0]) if G_el is not None else G[0], G_yield, color='#2a78d6', alpha=0.08, lw=0)
    if ep and y_el is not None:   # elastic first, the elasto-plastic (dashed) on top
        ax.plot(G_el, y_el, color='#2a78d6', lw=2, label=ylabel.split(',')[0] + ', elastic surface')
    ax.plot(G, y, color='k' if not ep else '#eb6834', lw=2, ls='-' if not ep else '--',
            label=main_label or (ylabel.split(',')[0] + sfx), zorder=3)
    for xe, ye, lab, kw in extra:   # components of the quantity
        ax.plot(xe, ye, label=lab, **kw)
    if ep and G_yield:
        ax.axvline(G_yield, color='k', lw=1, ls='-.', label=r'Yield (elastic range shaded), $\nabla T$ = %.0f K/m' % G_yield)
    if G_mean:
        ax.axvline(G_mean, color='0.4', lw=1, ls=':', label=r'Mean experimental $\overline{\nabla T}$ = %g K/m' % G_mean)
    txt = []
    if G[-1] >= G_mark:
        ax.axvline(G_mark, color='#b02a2a', lw=1, ls='--')
        v = np.interp(G_mark, G, y)
        ax.plot([G_mark], [v], 'o', color='#b02a2a', ms=6, zorder=5)
        txt.append(('elasto-plastic' if ep else 'elastic') + f' {v:.4g} {unit}')
        if ep and y_el is not None and G_el[-1] >= G_mark:
            ve = np.interp(G_mark, G_el, y_el)
            ax.plot([G_mark], [ve], 'o', color='#b02a2a', ms=6, zorder=5)
            txt.append(f'elastic {ve:.4g} {unit} ({100 * (v / ve - 1):+.1f} %)')
    if txt:
        ax.text(0.03, 0.97 if log_y else 0.62, r'$\nabla T = 10^5$ K/m:' + '\n' + '\n'.join(txt), transform=ax.transAxes,
                va='top', fontsize=9, bbox=dict(boxstyle='round', fc='white', ec='0.6'))
    ax.set_xscale('log')
    if log_y:
        ax.set_yscale('log')
    ax.set_xlabel(r'Thermal Gradient, $\mathrm{\nabla T \,\, [K.m^{-1}]}$', {'color': 'black', 'fontsize': 15})
    ax.set_ylabel(ylabel, {'color': 'black', 'fontsize': 15})
    plt.setp(ax.get_xticklabels(), fontsize=13)
    plt.setp(ax.get_yticklabels(), fontsize=13)
    ax.grid(True, which='both', linestyle='--', color='0.85', linewidth=0.5)
    ax.legend(fontsize=8.5, frameon=True, edgecolor='k', shadow=True, loc='lower right' if log_y else 'lower left')
    return fig, ax


def texts_inside_all():
    """Never prevents the figures from being shown: a figure that cannot be adjusted is left as it is."""
    try:
        import mpl_textguard   # noqa: F401  (off-screen texts with matplotlib >= 3.10, if the helper is present)
    except Exception:
        pass
    for n in plt.get_fignums():
        try:
            texts_inside(plt.figure(n))
        except Exception as e:
            print(f'figure_layout: figure {n} left unchanged ({type(e).__name__}: {e})')
