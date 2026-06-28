"""Shared postprocessing functions for coalescence simulations."""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

# Publication-quality matplotlib configuration
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['Computer Modern Roman']
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
matplotlib.rcParams['figure.dpi'] = 150
matplotlib.rcParams['lines.linewidth'] = 2.5

plt_settings = {
    'LabelFont': 28,
    'AxesFont': 22,
    'TitleFont': 28,
    'LegendFont': 18,
    'ColorbarFont': 22,
}


def style_axis(ax, xlabel=None, ylabel=None, title=None):
    """Apply publication-quality styling to an axis."""
    ax.tick_params(axis='both', which='major', labelsize=plt_settings['AxesFont'],
                   width=2, length=8, direction='out', pad=8)
    ax.tick_params(which='minor', width=1.5, length=4, direction='out')
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    ax.minorticks_on()
    ax.grid(True, alpha=0.3, linewidth=1)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=plt_settings['LabelFont'], labelpad=10)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=plt_settings['LabelFont'], labelpad=10)
    if title:
        ax.set_title(title, fontsize=plt_settings['TitleFont'], pad=15)


def find_neck(x, h):
    """Find coalescence neck (local minimum closest to x=0).

    Uses scipy.signal.find_peaks with prominence filtering to robustly
    identify local minima, then selects the one closest to x=0.

    Args:
        x: position array
        h: height array

    Returns:
        x0: neck position
        h0: neck height
    """
    sort_idx = np.argsort(x)
    x_sorted, h_sorted = x[sort_idx], h[sort_idx]

    # Find local minima (peaks in -h) with prominence filtering
    prominence = 0.01 * (h_sorted.max() - h_sorted.min())
    min_indices, _ = find_peaks(-h_sorted, prominence=prominence)

    if len(min_indices) == 0:
        # Fallback: interpolate at x=0
        return 0.0, np.interp(0.0, x_sorted, h_sorted)

    # Select minimum closest to x=0
    min_positions = x_sorted[min_indices]
    closest_idx = min_indices[np.argmin(np.abs(min_positions))]

    return x_sorted[closest_idx], h_sorted[closest_idx]


def find_drop_edge(x, h, hp, side='right'):
    """Find drop edge position where h drops below threshold.

    Args:
        x: position array
        h: height array
        hp: precursor film thickness
        side: 'right' or 'left'

    Returns:
        xe: drop edge position
    """
    sort_idx = np.argsort(x)
    x_sorted, h_sorted = x[sort_idx], h[sort_idx]

    threshold = max(2 * hp, 0.01)
    if side == 'right':
        mask = x_sorted > 0
    else:
        mask = x_sorted < 0

    x_side, h_side = x_sorted[mask], h_sorted[mask]
    below = h_side < threshold

    if np.any(below):
        return x_side[np.argmax(below)]
    return x_side[-1] if side == 'right' else x_side[0]
