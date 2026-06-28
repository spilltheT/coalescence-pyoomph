#!/usr/bin/env python3
"""
HvsT.py - Plot h(t) for droplet coalescence with surfactants.

Generates incremental publication-quality log-log plots for presentations:
- _v0: Empty axes (box only)
- _v1: Single case (β=0.8, Pe=1, θ=10°)
- _v2: All θ=10° cases (4 blue datasets)
- _v3: All data (8 datasets)

Visual encoding:
- Color: θ (Blue = 10°, Orange = 20°)
- Marker: β (Circle = 0.8, Square = 0.5)
- Shade: Pe (Light = 1, Dark = 10)

Converted from HvsT.m
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# =============================================================================
# Publication-quality plot settings (LaTeX fonts, thick spines)
# =============================================================================
plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'serif',
    'font.serif': ['Computer Modern Roman'],
    'font.size': 10,
    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'legend.fontsize': 9,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'axes.linewidth': 1.5,
    'xtick.major.width': 1.2,
    'ytick.major.width': 1.2,
    'xtick.major.size': 5,
    'ytick.major.size': 5,
    'xtick.minor.size': 3,
    'ytick.minor.size': 3,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'legend.frameon': False,
    'figure.dpi': 150,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})

# =============================================================================
# Output directories
# =============================================================================
root_dir = Path(__file__).parent
plots_dir = root_dir / 'HvsT_plots'
plots_dir.mkdir(exist_ok=True)

# =============================================================================
# Data loading
# =============================================================================
data_dir = root_dir / 'h_f_vals'

# Load all 8 datasets
min_h = []
for k in range(1, 9):
    file_path = data_dir / f'min_h_values_{k}.txt'
    df = pd.read_csv(file_path, delimiter='\t')
    min_h.append(df['Min_h'].values)

# =============================================================================
# Parameters
# =============================================================================
# Time arrays
t1 = np.linspace(0, 1000, 10001)  # For datasets 1-4 (theta = 10 deg)
t2 = np.linspace(0, 250, 2501)    # For datasets 5-8 (theta = 20 deg)

# t_0 correction for t1
t1 = t1 + 0.39

# Physical parameters for each dataset
beta = np.array([0.8, 0.8, 0.5, 0.5, 0.8, 0.8, 0.5, 0.5])
Pe = np.array([1, 10, 1, 10, 1, 10, 1, 10])
theta1 = 10 * np.pi / 180
theta2 = 20 * np.pi / 180
theta = np.array([theta1, theta1, theta1, theta1, theta2, theta2, theta2, theta2])
theta_deg = np.array([10, 10, 10, 10, 20, 20, 20, 20])

# Theoretical scaling for compensation
Gamma_0 = 0.8
h_a = (1 - beta * Gamma_0 / 2) * theta**4

# =============================================================================
# Visual encoding system
# =============================================================================
# Color by theta: Blue for 10°, Orange for 20°
# Shade by Pe: Light for Pe=1, Dark for Pe=10
colors_theta10 = {1: '#6baed6', 10: '#08519c'}   # Light/dark blue
colors_theta20 = {1: '#fd8d3c', 10: '#d94701'}   # Light/dark orange

# Marker by beta: Circle for 0.8, Square for 0.5
markers_beta = {0.8: 'o', 0.5: 's'}

def get_style(beta_val, pe_val, theta_val):
    """Return (color, marker) based on parameters."""
    if theta_val == 10:
        color = colors_theta10[pe_val]
    else:
        color = colors_theta20[pe_val]
    marker = markers_beta[beta_val]
    return color, marker

# =============================================================================
# Plot configuration: which datasets to plot and their ranges
# =============================================================================
# Format: (dataset_index, time_array, end_index, beta, Pe, theta_deg)
plot_config = [
    (0, t1, 10001, 0.8, 1, 10),
    (1, t1, 10001, 0.8, 10, 10),
    (2, t1, 10001, 0.5, 1, 10),
    (3, t1, 10001, 0.5, 10, 10),
    (4, t2, 2501, 0.8, 1, 20),
    (5, t2, 2501, 0.8, 10, 20),
    (6, t2, 2501, 0.5, 1, 20),
    (7, t2, 2501, 0.5, 10, 20),
]

# Marker settings
marker_size = 8

def log_subsample(t_arr, y_arr, n_points=200):
    """Subsample data uniformly in log-space for even distribution on log-log plots."""
    # Find valid range (t > 0)
    valid = t_arr > 0
    t_valid = t_arr[valid]
    y_valid = y_arr[valid]

    if len(t_valid) == 0:
        return t_arr, y_arr

    # Generate logarithmically spaced indices
    log_t_min, log_t_max = np.log10(t_valid.min()), np.log10(t_valid.max())
    log_t_targets = np.linspace(log_t_min, log_t_max, n_points)

    # Find nearest index for each target
    indices = []
    for log_t in log_t_targets:
        idx = np.argmin(np.abs(np.log10(t_valid) - log_t))
        if idx not in indices:  # Avoid duplicates
            indices.append(idx)

    indices = np.array(indices)
    return t_valid[indices], y_valid[indices]

# =============================================================================
# Incremental build stages
# =============================================================================
stages = [
    ('v0', []),                       # Empty box
    ('v1', [0]),                      # β=0.8, Pe=1, θ=10°
    ('v2', [0, 1, 2, 3]),             # All θ=10° (blues)
    ('v3', [0, 1, 2, 3, 4, 5, 6, 7]), # All data
]

# =============================================================================
# Plotting functions
# =============================================================================
def make_uncompensated_plot(indices_to_plot, show_refline=True):
    """Create uncompensated h_min(t) plot with specified datasets."""
    fig, ax = plt.subplots(figsize=(5, 4))

    # Set log scale explicitly (needed for empty plots)
    ax.set_xscale('log')
    ax.set_yscale('log')

    for idx, t_arr, end_idx, b, p, th in plot_config:
        if idx in indices_to_plot:
            color, marker = get_style(b, p, th)
            label = rf'$\beta={b}$, $\mathrm{{Pe}}={p}$, $\theta={th}^\circ$'
            # Log-subsample for even distribution on log-log plot
            t_plot, h_plot = log_subsample(t_arr[1:end_idx], min_h[idx][1:end_idx])
            ax.loglog(
                t_plot, h_plot,
                marker, markersize=marker_size,
                markerfacecolor=color, markeredgecolor='none',
                alpha=0.8, label=label
            )

    # Reference power law: h ~ t
    if show_refline and len(indices_to_plot) > 0:
        t_ref = np.logspace(np.log10(0.1), np.log10(500), 100)
        ax.loglog(t_ref, 0.000872 * t_ref**1, 'k-', linewidth=2.5, label=r'$\sim t$')

    ax.set_xlabel(r'$t$', fontsize=14)
    ax.set_ylabel(r'$h_0$', fontsize=14)

    # Set consistent axis limits for all versions
    ax.set_xlim([0.1, 1000])
    ax.set_ylim([1e-4, 0.2])

    if len(indices_to_plot) > 0:
        ax.legend(loc='lower right', fontsize=7, handletextpad=0.3, labelspacing=0.3)

    fig.tight_layout()
    return fig, ax


def make_compensated_plot(indices_to_plot, show_refline=True):
    """Create compensated h_min(t)/h_a plot with specified datasets."""
    fig, ax = plt.subplots(figsize=(5, 4))

    # Set log scale explicitly (needed for empty plots)
    ax.set_xscale('log')
    ax.set_yscale('log')

    for idx, t_arr, end_idx, b, p, th in plot_config:
        if idx in indices_to_plot:
            color, marker = get_style(b, p, th)
            label = rf'$\beta={b}$, $\mathrm{{Pe}}={p}$, $\theta={th}^\circ$'
            # Log-subsample for even distribution on log-log plot
            t_plot, h_plot = log_subsample(t_arr[1:end_idx], min_h[idx][1:end_idx] / h_a[idx])
            ax.loglog(
                t_plot, h_plot,
                marker, markersize=marker_size,
                markerfacecolor=color, markeredgecolor='none',
                alpha=0.8, label=label
            )

    # Reference power law: h/h_a ~ t
    if show_refline and len(indices_to_plot) > 0:
        t_ref = np.logspace(np.log10(0.1), np.log10(500), 100)
        ax.loglog(t_ref, 0.272 * t_ref**1, 'k-', linewidth=2.5, label=r'$\sim t$')

    ax.set_xlabel(r'$t$', fontsize=14)
    ax.set_ylabel(r'$h_0 / \left(1-\frac{\beta\Gamma_0}{2}\right)\theta^4$', fontsize=14)

    # Set consistent axis limits for all versions
    ax.set_xlim([0.1, 1000])
    ax.set_ylim([1e-2, 1e2])

    if len(indices_to_plot) > 0:
        ax.legend(loc='lower right', fontsize=7, handletextpad=0.3, labelspacing=0.3)

    fig.tight_layout()
    return fig, ax


# =============================================================================
# Generate all plots
# =============================================================================
for stage_name, indices in stages:
    # Uncompensated (no reference line in staged versions)
    fig, ax = make_uncompensated_plot(indices, show_refline=False)
    fig.savefig(plots_dir / f'HvsT_uncompensated_{stage_name}.pdf')
    plt.close(fig)
    print(f"Saved: HvsT_plots/HvsT_uncompensated_{stage_name}.pdf")

    # Compensated (no reference line in staged versions)
    fig, ax = make_compensated_plot(indices, show_refline=False)
    fig.savefig(plots_dir / f'HvsT_compensated_{stage_name}.pdf')
    plt.close(fig)
    print(f"Saved: HvsT_plots/HvsT_compensated_{stage_name}.pdf")

# =============================================================================
# Final outputs in root directory
# =============================================================================
# Uncompensated (all data)
fig, ax = make_uncompensated_plot([0, 1, 2, 3, 4, 5, 6, 7])
fig.savefig(root_dir / 'HvsT_uncompensated.pdf')
plt.close(fig)
print("Saved: HvsT_uncompensated.pdf")

# Compensated (all data)
fig, ax = make_compensated_plot([0, 1, 2, 3, 4, 5, 6, 7])
fig.savefig(root_dir / 'HvsT_compensated.pdf')
plt.close(fig)
print("Saved: HvsT_compensated.pdf")
