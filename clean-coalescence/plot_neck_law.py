#!/usr/bin/env python3
"""Neck-height law for clean sessile-drop coalescence.

In the viscous lubrication limit the bridge neck grows as

    h0(t) = 0.272 theta^4 t

(Hernandez-Sanchez et al. 2012). Compensating the neck history of every reduced
run in ``data/`` by 0.272 theta^4 collapses all contact angles onto the single
universal line h0/(0.272 theta^4) = t. Each angle peels off only once the bridge
feels the finite drop.

    python plot_neck_law.py        # writes neck_law.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
V_STAR = 0.272
HP = 1e-4         # precursor film thickness used in the simulations

COLORS = {1: "#D81B7A", 2: "#7D3C98", 3: "#1B9E77",
          5: "#C0392B", 10: "#1F6FA8", 20: "#E08214"}
MARKS = {1: "P", 2: "v", 3: "D", 5: "^", 10: "s", 20: "o"}


def _configure(use_tex=True):
    matplotlib.rcParams.update({
        "font.family": "serif", "font.serif": ["Computer Modern Roman"],
        "font.size": 11, "axes.linewidth": 1.0, "text.usetex": use_tex, "axes.unicode_minus": False,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
    })
    if not use_tex:
        matplotlib.rcParams["mathtext.fontset"] = "cm"


def _subsample_log(t, y, n=24):
    m = t > 0
    t, y = t[m], y[m]
    g = np.logspace(np.log10(t.min()), np.log10(t.max()), n)
    idx = np.unique([int(np.argmin(np.abs(t - q))) for q in g])
    return t[idx], y[idx]


def build(use_tex=True):
    _configure(use_tex)
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    ax.set_xscale("log")
    ax.set_yscale("log")
    bundles = sorted(DATA.glob("*.npz"),
                     key=lambda p: float(np.load(p, allow_pickle=True)["theta_deg"]),
                     reverse=True)
    handles = []
    tmax = 1.0
    for p in bundles:
        d = np.load(p, allow_pickle=True)
        th = int(round(float(d["theta_deg"])))
        thr = np.deg2rad(th)
        color, mk = COLORS.get(th, "#444"), MARKS.get(th, "o")
        t, h = d["h0t"][:, 0], d["h0t"][:, 1]
        keep = h > 1.5 * HP
        if not keep.any():
            print(f"skipping {p.name}: no neck heights above {1.5 * HP:g} "
                  f"(truncated/too-short run)")
            continue
        t, y = t[keep], h[keep] / (V_STAR * thr**4)
        ts, ys = _subsample_log(t, y)
        ax.plot(ts, ys, mk, ms=5.4, mfc=color, mec="white", mew=0.3,
                ls="none", zorder=4)
        tmax = max(tmax, t.max())
        handles.append(plt.Line2D([], [], color=color, marker=mk, ls="none",
                                  ms=5.4, mec="white", mew=0.3,
                                  label=rf"$\theta={th}^\circ$"))
    tref = np.logspace(np.log10(0.2), np.log10(tmax * 1.3), 100)
    ax.plot(tref, tref, "k-", lw=1.7, zorder=6)
    handles.append(plt.Line2D([], [], color="k", lw=1.7,
                              label=r"$h_0 = 0.272\,\theta^4 t$"))
    ax.set_xlim(0.2, tmax * 2)
    ax.set_ylim(0.02, tmax * 2)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$h_0 \,/\, (0.272\,\theta^4)$")
    ax.legend(handles=handles, loc="lower right", ncol=2, frameon=True,
              fontsize=8, handletextpad=0.3, columnspacing=0.8)
    fig.tight_layout()
    return fig


def main():
    try:
        fig = build(use_tex=True)
    except RuntimeError:
        plt.close("all")
        fig = build(use_tex=False)
    out = HERE / "neck_law.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
