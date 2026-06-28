#!/usr/bin/env python3
"""
Presentation-quality single-panel figure for Marangoni cancellation analysis.

Shows chi_max vs Pe for different beta values using RdYlBu_r colormap.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os

# LaTeX configuration for publication-quality output
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['Computer Modern Roman']
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

# Font sizes (all set to 9)
plt_settings = {
    'LabelFont': 9,
    'AxesFont': 9,
    'TitleFont': 9,
    'LegendFont': 9,
    'ColorbarFont': 9,
}

# Figure dimensions for single panel
MM_TO_INCH = 1 / 25.4
PANEL_WIDTH_MM = 70
PANEL_HEIGHT_MM = 55


def style_axis(ax):
    """Apply publication-quality styling to an axis."""
    ax.tick_params(axis='both', which='major', labelsize=plt_settings['AxesFont'],
                   width=0.8, length=4, direction='in', pad=3)
    ax.tick_params(which='minor', width=0.5, length=2, direction='in')
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    ax.minorticks_on()
    ax.grid(True, alpha=0.25, linewidth=0.4)


def load_csv_data(filepath):
    """Load chi_max data from CSV file."""
    data = {}
    with open(filepath, 'r') as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split(',')
            beta = float(parts[0])
            pe = float(parts[1])
            chi_max = float(parts[2]) if parts[2] not in ['0', '1E-04', '0.0001', '0.0002'] else float(parts[2])

            if beta not in data:
                data[beta] = {'Pe': [], 'chiMax': []}
            data[beta]['Pe'].append(pe)
            data[beta]['chiMax'].append(chi_max)

    # Convert to numpy arrays
    for beta in data:
        data[beta]['Pe'] = np.array(data[beta]['Pe'])
        data[beta]['chiMax'] = np.array(data[beta]['chiMax'])

    return data


def main():
    # Data path
    csv_file = 'marangoni_cancellation_beta_Pe_chiMin.csv'

    # Load data
    print("Loading chi_max data...")
    csv_data = load_csv_data(csv_file)

    # Create single-panel figure
    fig_width = PANEL_WIDTH_MM * MM_TO_INCH
    fig_height = PANEL_HEIGHT_MM * MM_TO_INCH
    fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))
    fig.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95)
    fig.set_facecolor('white')

    # =========================================================================
    # chi_max vs Pe for different beta values
    # =========================================================================
    print("Plotting chi_max vs Pe...")

    # Marker styles for different beta values
    markers = ['o', 's', 'D', '^', 'v', 'p']
    # Selected beta values (subset for clarity)
    beta_values = [0.05, 0.1, 0.15, 0.2, 0.25, 0.5]

    # Setup colormap for beta using RdYlBu_r (blue->yellow->red for low->high beta)
    beta_norm = plt.Normalize(min(beta_values), max(beta_values))
    beta_cmap = cm.RdYlBu_r

    for i, beta in enumerate(beta_values):
        if beta not in csv_data:
            print(f"Warning: beta = {beta} not found in data, skipping")
            continue
        d = csv_data[beta]
        color = beta_cmap(beta_norm(beta))
        marker = markers[i % len(markers)]

        ax.semilogx(d['Pe'], d['chiMax'],
                    marker=marker,
                    linestyle='None',
                    markersize=5,
                    markerfacecolor=color,
                    markeredgecolor='k',
                    markeredgewidth=0.4,
                    label=rf'$\beta = {beta}$',
                    zorder=3)

    # Axis configuration
    ax.set_xlim(20, 1200)
    ax.set_ylim(0, 0.35)
    ax.set_xlabel(r'$Pe$', fontsize=plt_settings['LabelFont'], labelpad=3)
    ax.set_ylabel(r'$\chi_\mathrm{max}$', fontsize=plt_settings['LabelFont'], labelpad=3)

    # Legend
    ax.legend(loc='upper left', frameon=True, framealpha=0.9,
              fontsize=plt_settings['LegendFont'], ncol=1,
              handlelength=0.8, markerscale=0.8, handletextpad=0.3,
              borderpad=0.3, labelspacing=0.2)

    # Style axis
    style_axis(ax)

    # =========================================================================
    # Final adjustments and save
    # =========================================================================
    plt.tight_layout()

    output_file = 'fig_chi_analysis_PresentationOnly.pdf'
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()

    print(f"\nFigure saved to: {output_file}")


if __name__ == '__main__':
    main()
