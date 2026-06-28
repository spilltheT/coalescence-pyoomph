#!/usr/bin/env python3
"""
Publication-quality two-panel figure for Marangoni cancellation analysis.

Panel (a): chi(t) vs time for different Pe values (beta=0.1)
Panel (b): chi_max vs Pe for all beta values
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

# Figure dimensions: 135mm total width, 1:1 panel aspect ratio
MM_TO_INCH = 1 / 25.4
FIGURE_WIDTH_MM = 135
PANEL_HEIGHT_MM = 50  # Square panels


def style_axis(ax):
    """Apply publication-quality styling to an axis."""
    ax.tick_params(axis='both', which='major', labelsize=plt_settings['AxesFont'],
                   width=0.8, length=4, direction='in', pad=3)
    ax.tick_params(which='minor', width=0.5, length=2, direction='in')
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    ax.minorticks_on()
    ax.grid(True, alpha=0.25, linewidth=0.4)


def load_chi_data(folder, pe_values):
    """Load chi(t) data for all Pe values."""
    all_data = {}
    for pe in pe_values:
        filepath = os.path.join(folder, f'chi_Pe_{pe}.txt')
        if not os.path.exists(filepath):
            print(f"Warning: {filepath} not found, skipping")
            continue

        data = np.loadtxt(filepath)
        all_data[pe] = {
            'time': data[:, 0],
            'chi': data[:, 1],
        }
    return all_data


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
    # Pe values for Panel 1 (start at 250 where chi becomes non-zero)
    pe_values = [250, 300, 350, 400, 500, 600, 700, 800, 900, 1000]

    # Data paths
    chi_folder = 'chi-vals-beta0.1'
    csv_file = 'marangoni_cancellation_beta_Pe_chiMin.csv'

    # Load data
    print("Loading chi(t) data...")
    chi_data = load_chi_data(chi_folder, pe_values)

    print("Loading chi_max data...")
    csv_data = load_csv_data(csv_file)

    # Create figure with 2 side-by-side panels (135mm total width)
    fig_width = FIGURE_WIDTH_MM * MM_TO_INCH
    fig_height = PANEL_HEIGHT_MM * MM_TO_INCH + 0.4  # Extra space for labels
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fig_width, fig_height))
    fig.subplots_adjust(left=0.08, right=0.92, bottom=0.15, top=0.95, wspace=0.35)
    fig.set_facecolor('white')

    # =========================================================================
    # Panel 1: chi(t) vs time with Pe colorbar
    # =========================================================================
    print("Plotting Panel (a): chi(t) vs time...")

    # Setup colormap for Pe (log scale)
    log_pe = np.log10(pe_values)
    norm = plt.Normalize(log_pe.min(), log_pe.max())
    cmap = cm.viridis

    for pe in pe_values:
        if pe in chi_data:
            d = chi_data[pe]
            color = cmap(norm(np.log10(pe)))
            ax1.plot(d['time'], d['chi'], color=color, linewidth=0.8)

    # Axis configuration
    ax1.set_ylim(0, 0.2)
    ax1.set_xlim(0, 100)
    ax1.set_xlabel(r'$t$', fontsize=plt_settings['LabelFont'], labelpad=3)
    ax1.set_ylabel(r'$\chi$', fontsize=plt_settings['LabelFont'], labelpad=3)

    # Colorbar for Pe
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax1, pad=0.02, aspect=15)
    cbar.set_label(r'$\log_{10}(Pe)$', fontsize=plt_settings['ColorbarFont'], labelpad=3)
    cbar.ax.tick_params(labelsize=plt_settings['AxesFont'], width=0.5, length=2)

    # Style axis
    style_axis(ax1)

    # Add chi_max indicator for Pe = 1000
    if 1000 in chi_data:
        chi_max_pe1000 = np.max(chi_data[1000]['chi'])
        # Draw short horizontal line (black for neutral contrast with viridis)
        ax1.hlines(chi_max_pe1000, 30, 70, colors='k', linewidth=1.2, zorder=5)
        # Add label
        ax1.text(72, chi_max_pe1000, r'$\chi_\mathrm{max}$', fontsize=plt_settings['LabelFont'],
                 va='center', ha='left')

    # Panel label (outside, top left)
    ax1.text(-0.25, 1.05, r'(a)', transform=ax1.transAxes, fontsize=plt_settings['TitleFont'],
             verticalalignment='bottom', fontweight='bold')

    # =========================================================================
    # Panel 2: chi_max vs Pe for different beta values
    # =========================================================================
    print("Plotting Panel (b): chi_max vs Pe...")

    # Marker styles for different beta values
    markers = ['o', 's', 'D', '^', 'v']
    # Only use selected beta values
    beta_values = [0.05, 0.1, 0.15, 0.25, 0.5]

    # Setup colormap for beta (plasma: perceptually uniform, all colors visible)
    beta_norm = plt.Normalize(min(beta_values), max(beta_values))
    beta_cmap = cm.plasma

    for i, beta in enumerate(beta_values):
        if beta not in csv_data:
            print(f"Warning: beta = {beta} not found in data, skipping")
            continue
        d = csv_data[beta]
        color = beta_cmap(beta_norm(beta))
        marker = markers[i % len(markers)]

        ax2.semilogx(d['Pe'], d['chiMax'],
                     marker=marker,
                     linestyle='None',
                     markersize=5,
                     markerfacecolor=color,
                     markeredgecolor='k',
                     markeredgewidth=0.4,
                     label=rf'$\beta = {beta}$',
                     zorder=3)

    # Axis configuration
    ax2.set_xlim(20, 1200)
    ax2.set_ylim(0, 0.35)
    ax2.set_xlabel(r'$Pe$', fontsize=plt_settings['LabelFont'], labelpad=3)
    ax2.set_ylabel(r'$\chi_\mathrm{max}$', fontsize=plt_settings['LabelFont'], labelpad=3)

    # Legend
    ax2.legend(loc='upper left', frameon=True, framealpha=0.9,
               fontsize=plt_settings['LegendFont'], ncol=1,
               handlelength=0.8, markerscale=0.8, handletextpad=0.3,
               borderpad=0.3, labelspacing=0.2)

    # Style axis
    style_axis(ax2)

    # Panel label (outside, top left)
    ax2.text(-0.15, 1.05, r'(b)', transform=ax2.transAxes, fontsize=plt_settings['TitleFont'],
             verticalalignment='bottom', fontweight='bold')

    # =========================================================================
    # Final adjustments and save
    # =========================================================================
    plt.tight_layout()

    output_file = 'fig_chi_analysis.pdf'
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()

    print(f"\nFigure saved to: {output_file}")


if __name__ == '__main__':
    main()
