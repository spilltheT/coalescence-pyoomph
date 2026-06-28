#!/usr/bin/env python3
"""Self-similar collapse of clean sessile-drop coalescence.

For each reduced run in ``data/`` this plots the rescaled bridge profiles

    H = h / h0(t)      versus      xi = theta x / h0(t)

(measured-neck-height variables of Hernandez-Sanchez et al. 2012), overlaid
with the self-similar master curve computed by ``similarity_solution.py``.
Profiles taken at different times and different contact angles fall onto the
single master curve: the early bridge is self-similar and forgets its initial
geometry.

    python plot_collapse.py        # writes self_similar_collapse.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
MASTER = HERE / "similarity_master.csv"
XI = 5.0          # plot range in the similarity variable
NCURVES = 8       # snapshots overlaid per run

COLORS = {1: "#D81B7A", 2: "#7D3C98", 3: "#1B9E77",
          5: "#C0392B", 10: "#1F6FA8", 20: "#E08214"}


def _configure(use_tex=True):
    matplotlib.rcParams.update({
        "font.family": "serif", "font.serif": ["Computer Modern Roman"],
        "font.size": 11, "axes.linewidth": 1.0, "text.usetex": use_tex, "axes.unicode_minus": False,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
    })
    if not use_tex:
        matplotlib.rcParams["mathtext.fontset"] = "cm"


def build(use_tex=True):
    _configure(use_tex)
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    bundles = sorted(DATA.glob("*.npz"),
                     key=lambda p: float(np.load(p, allow_pickle=True)["theta_deg"]),
                     reverse=True)
    handles = []
    for k, p in enumerate(bundles):
        d = np.load(p, allow_pickle=True)
        th = int(round(float(d["theta_deg"])))
        color = COLORS.get(th, "#444444")
        n = len(d["coll_t"])
        idx = np.unique(np.linspace(0, n - 1, min(n, NCURVES)).astype(int))
        lw, ls = (1.7, "-") if k == 0 else (1.3, (0, (5, 2)))
        for i in idx:
            ax.plot(d["coll_xi"][i], d["coll_H"][i], color=color, lw=lw, ls=ls,
                    alpha=0.95, solid_capstyle="round", zorder=3 + k)
        handles.append(plt.Line2D([], [], color=color, lw=2.2, ls=ls,
                                  label=rf"$\theta={th}^\circ$"))
    if MASTER.exists():
        m = np.loadtxt(MASTER, delimiter=",")
        sel = np.abs(m[:, 0]) <= XI
        ax.plot(m[sel, 0], m[sel, 1], color="k", lw=1.8, zorder=9)
        handles.append(plt.Line2D([], [], color="k", lw=1.8, label="theory"))
    ax.set_xlim(-XI, XI)
    ax.set_ylim(0, XI)
    ax.set_xlabel(r"$\xi = \theta x / h_0(t)$")
    ax.set_ylabel(r"$\mathcal{H} = h / h_0(t)$")
    ax.legend(handles=handles, loc="upper center", ncol=len(handles),
              frameon=False, fontsize=9, handlelength=1.4, columnspacing=0.9)
    fig.tight_layout()
    return fig


def main():
    try:
        fig = build(use_tex=True)
    except RuntimeError:
        plt.close("all")
        fig = build(use_tex=False)
    out = HERE / "self_similar_collapse.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
