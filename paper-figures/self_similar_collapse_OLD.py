#!/usr/bin/env python3
"""
self_similar_collapse.py - Plot self-similar collapse of height profiles h/h_0(t) vs x/h_0(t).

Demonstrates that height profiles near the coalescence neck DO NOT collapse to a universal
self-similar shape when rescaled by h_0(t) alone (without theta correction).

Plots two datasets (θ=10° and θ=20°) to show the scaling fails without theta factor.

Scaling:
    y-axis: h / h_0(t)
    x-axis: x / h_0(t)  (OLD - missing theta factor)

Usage:
    python self_similar_collapse.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pathlib import Path
import re

# =============================================================================
# Publication-quality plot settings (LaTeX fonts, thick spines)
# =============================================================================
plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'serif',
    'font.serif': ['Computer Modern Roman'],
    'font.size': 10,
    'axes.labelsize': 14,
    'axes.titlesize': 12,
    'legend.fontsize': 9,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'axes.linewidth': 2.0,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'xtick.major.size': 6,
    'ytick.major.size': 6,
    'xtick.minor.size': 3,
    'ytick.minor.size': 3,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'legend.frameon': False,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})


def load_domain_file(filepath):
    """
    Load a domain snapshot file.

    Returns:
        time: float, simulation time
        x: array, spatial coordinates
        h: array, height profile
    """
    with open(filepath, 'r') as f:
        header = f.readline()

    # Extract time from header: "# coordinate_x	h	p	Gamma	@time=40.0"
    match = re.search(r'@time=([0-9.eE+-]+)', header)
    if match:
        time = float(match.group(1))
    else:
        time = 0.0

    # Load data (skip header)
    data = np.loadtxt(filepath, skiprows=1)
    x = data[:, 0]
    h = data[:, 1]

    return time, x, h


def find_neck_minimum(x, h):
    """
    Find the neck minimum (h_0) near x=0.

    Returns:
        h0: minimum height value
        x0: x-position of minimum
    """
    # Focus on region near x=0
    mask = np.abs(x) < 1.0
    if np.sum(mask) < 3:
        mask = np.abs(x) < 2.0

    x_near = x[mask]
    h_near = h[mask]

    idx_min = np.argmin(h_near)
    h0 = h_near[idx_min]
    x0 = x_near[idx_min]

    return h0, x0


def load_snapshots(data_dir, t_min, t_max):
    """Load all snapshots from a directory within time range."""
    data_path = Path(__file__).parent / data_dir
    if not data_path.exists():
        data_path = Path(data_dir)

    domain_files = sorted(data_path.glob('domain_*.txt'))

    if not domain_files:
        print(f"Error: No domain_*.txt files found in {data_path}")
        return []

    print(f"Found {len(domain_files)} domain files in {data_path}")

    snapshots = []
    for filepath in domain_files:
        time, x, h = load_domain_file(filepath)
        if t_min <= time <= t_max:
            h0, x0 = find_neck_minimum(x, h)
            snapshots.append({
                'time': time,
                'x': x,
                'h': h,
                'h0': h0,
                'x0': x0
            })

    snapshots.sort(key=lambda s: s['time'])
    print(f"  {len(snapshots)} snapshots in time range [{t_min}, {t_max}]")

    return snapshots


def main():
    # Configuration
    x_lim, h_lim = 20.0, 4.0
    output_file = 'self_similar_collapse.pdf'

    # Dataset configurations: (data_dir, theta_deg, colormap, t_min, t_max, label)
    datasets = [
        ('hx-Pe1-beta0p1-theta10', 10.0, cm.Blues, 2.0, 8.0, r'$\theta=10^\circ$'),
        ('hx-Pe1-beta0p1-theta20', 20.0, cm.hot_r, 0.5, 2.0, r'$\theta=20^\circ$'),
    ]

    # Create figure
    fig, ax = plt.subplots(figsize=(7, 5))

    # Store normalization for each dataset (for colorbars)
    norms = []

    # Plot each dataset
    for data_dir, theta_deg, cmap, t_min, t_max, label in datasets:
        print(f"\nProcessing {label}:")
        print(f"  Time range: {t_min} < t < {t_max}")

        norm = plt.Normalize(t_min, t_max)
        norms.append((norm, cmap, t_min, t_max, label))

        snapshots = load_snapshots(data_dir, t_min, t_max)
        if not snapshots:
            continue

        # Plot each snapshot
        for snap in snapshots:
            t = snap['time']
            x0 = snap['x0']
            h0 = snap['h0']

            # OLD scaling: just h_0, NO theta factor on x
            x_scaled = (snap['x'] - x0) / h0
            h_scaled = snap['h'] / h0

            color = cmap(norm(t))
            ax.plot(x_scaled, h_scaled, 'o',
                    markersize=2, color=color, alpha=0.7,
                    markeredgecolor='none')

    # Axis labels
    ax.set_xlabel(r'$\displaystyle\frac{x}{h_0(t)}$')
    ax.set_ylabel(r'$\displaystyle\frac{h(x,t)}{h_0(t)}$')

    # Axis limits
    ax.set_xlim(-x_lim, x_lim)
    ax.set_ylim(0, h_lim)

    # Create two colorbars side by side
    divider = make_axes_locatable(ax)

    # First colorbar (θ=10°, Blues)
    norm1, cmap1, t_min1, t_max1, label1 = norms[0]
    cax1 = divider.append_axes("right", size="4%", pad=0.1)
    sm1 = cm.ScalarMappable(cmap=cmap1, norm=norm1)
    sm1.set_array([])
    cbar1 = plt.colorbar(sm1, cax=cax1)
    cbar1.set_label(r'$t$ ($\theta=10^\circ$)', fontsize=11)
    cbar1.ax.tick_params(labelsize=10)

    # Second colorbar (θ=20°, hot_r)
    norm2, cmap2, t_min2, t_max2, label2 = norms[1]
    cax2 = divider.append_axes("right", size="4%", pad=0.4)
    sm2 = cm.ScalarMappable(cmap=cmap2, norm=norm2)
    sm2.set_array([])
    cbar2 = plt.colorbar(sm2, cax=cax2)
    cbar2.set_label(r'$t$ ($\theta=20^\circ$)', fontsize=11)
    cbar2.ax.tick_params(labelsize=10)

    fig.tight_layout()

    # Save
    output_path = Path(__file__).parent / output_file
    fig.savefig(output_path)
    print(f"\nSaved: {output_path}")

    plt.close(fig)


if __name__ == '__main__':
    main()
