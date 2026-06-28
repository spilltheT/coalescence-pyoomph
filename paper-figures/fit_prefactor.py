#!/usr/bin/env python3
"""
fit_prefactor.py - Fit the prefactor A for x_f = A * (t - t_0)^(3/2) scaling.

Automatically detects the valid scaling regime by computing local slopes
in log-log space and finding regions where slope ≈ 1.5. Uses case-dependent
t_0 (minimum of valid t range for each dataset) to shift the time origin.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

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
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})


def compute_local_slope(t, y, window=5):
    """
    Compute local logarithmic slope d[log(y)]/d[log(t)] using rolling window.

    Parameters
    ----------
    t : array-like
        Time values (must be positive)
    y : array-like
        Data values (must be positive)
    window : int
        Window size for rolling linear regression

    Returns
    -------
    slopes : ndarray
        Local slopes at each point (NaN at edges)
    """
    log_t = np.log(t)
    log_y = np.log(y)
    n = len(t)
    slopes = np.full(n, np.nan)

    half_w = window // 2
    for i in range(half_w, n - half_w):
        idx = slice(i - half_w, i + half_w + 1)
        # Linear regression: log(y) = slope * log(t) + intercept
        slope, _, _, _, _ = stats.linregress(log_t[idx], log_y[idx])
        slopes[i] = slope

    return slopes


def find_valid_indices(slopes, target=1.5, tol=0.5, min_points=10):
    """
    Find indices where slope is within tolerance of target.

    Parameters
    ----------
    slopes : array-like
        Local slope values
    target : float
        Target slope (default 1.5 for t^(3/2))
    tol : float
        Tolerance for valid slope
    min_points : int
        Minimum number of valid points required

    Returns
    -------
    valid_mask : ndarray
        Boolean mask for valid indices
    """
    valid = np.abs(slopes - target) < tol
    valid[np.isnan(slopes)] = False

    # Only keep if we have enough points
    if np.sum(valid) < min_points:
        return np.zeros_like(valid, dtype=bool)

    return valid


# =============================================================================
# Data loading (from XvsT.py)
# =============================================================================
data_dir = Path(__file__).parent / 'h_f_vals'

# Load all 8 datasets
min_x = []
for k in range(1, 9):
    file_path = data_dir / f'min_h_values_{k}.txt'
    df = pd.read_csv(file_path, delimiter='\t')
    min_x.append(df['Min_x'].values)

# =============================================================================
# Parameters (from XvsT.py)
# =============================================================================
# Time arrays
t1 = np.linspace(0, 1000, 10001)  # For datasets 1-4 (theta = 10 deg)
t2 = np.linspace(0, 250, 2501)    # For datasets 5-8 (theta = 20 deg)

# Physical parameters for each dataset
beta = np.array([0.8, 0.8, 0.5, 0.5, 0.8, 0.8, 0.5, 0.5])
Pe = np.array([1, 10, 1, 10, 1, 10, 1, 10])
theta1 = 10 * np.pi / 180
theta2 = 20 * np.pi / 180
theta = np.array([theta1, theta1, theta1, theta1, theta2, theta2, theta2, theta2])

# Theoretical scaling for compensation
Gamma_0 = 0.8
c_v = 0.272 * (1 - beta * Gamma_0 / 2) * theta**4
x_a = (c_v * beta * Gamma_0 * np.sqrt(Pe))
x_a = x_a / (3 * np.sqrt(np.pi))

# =============================================================================
# Colors (Tableau-10)
# =============================================================================
colors = [
    '#1f77b4',  # blue
    '#ff7f0e',  # orange
    '#2ca02c',  # green
    '#d62728',  # red
    '#9467bd',  # purple
    '#8c564b',  # brown
    '#e377c2',  # pink
    '#7f7f7f',  # gray
]

# =============================================================================
# Plot configuration (from XvsT.py)
# =============================================================================
plot_config = [
    (0, t1, 3000, r'$\beta=0.8$, $\mathrm{Pe}=1$, $\theta=10^\circ$'),
    (1, t1, 3000, r'$\beta=0.8$, $\mathrm{Pe}=10$, $\theta=10^\circ$'),
    (2, t1, 3000, r'$\beta=0.5$, $\mathrm{Pe}=1$, $\theta=10^\circ$'),
    (3, t1, 3000, r'$\beta=0.5$, $\mathrm{Pe}=10$, $\theta=10^\circ$'),
    (4, t2, 420, r'$\beta=0.8$, $\mathrm{Pe}=1$, $\theta=20^\circ$'),
    (5, t2, 550, r'$\beta=0.8$, $\mathrm{Pe}=10$, $\theta=20^\circ$'),
    (6, t2, 300, r'$\beta=0.5$, $\mathrm{Pe}=1$, $\theta=20^\circ$'),
    (7, t2, 450, r'$\beta=0.5$, $\mathrm{Pe}=10$, $\theta=20^\circ$'),
]

# Marker settings
marker_size = 8

# =============================================================================
# Global fit: collect valid points from all datasets
# =============================================================================
print("=" * 60)
print("Fitting prefactor for x_f = A * (t - t_0)^(3/2) scaling")
print("=" * 60)

all_valid_t_shifted = []  # Will store (t - t₀) values
all_valid_xf = []
dataset_t0 = {}  # Store t₀ for each dataset
dataset_t_range = {}  # Store (t_min, t_max) for valid range

for idx, t_arr, end_idx, label in plot_config:
    # Get compensated data
    t_data = t_arr[1:end_idx]
    xf_data = min_x[idx][1:end_idx] / x_a[idx]

    # Filter out any non-positive values
    mask_pos = (t_data > 0) & (xf_data > 0)
    t_data = t_data[mask_pos]
    xf_data = xf_data[mask_pos]

    # Compute local slopes
    slopes = compute_local_slope(t_data, xf_data, window=7)

    # Find valid region (slope ≈ 1.5)
    valid = find_valid_indices(slopes, target=1.5, tol=0.5, min_points=10)

    n_valid = np.sum(valid)
    print(f"\nDataset {idx}: {label}")
    print(f"  Valid points: {n_valid} / {len(t_data)}")

    if n_valid > 0:
        t_valid = t_data[valid]
        xf_valid = xf_data[valid]

        # Case-dependent t₀: minimum of valid t range
        t0 = t_valid.min()
        dataset_t0[idx] = t0
        if t_valid.max() < 50.0:
            temp = t_valid.max()
        else:
            temp = 50.0
        dataset_t_range[idx] = (t_valid.min(), temp)

        mask_tmax = t_valid <= temp
        t_valid = t_valid[mask_tmax]
        xf_valid = xf_valid[mask_tmax]

        print(f"  t_0 = {t0:.2f}")
        print(f"  Valid t range: [{t_valid.min():.2f}, {t_valid.max():.2f}]")
        print(f"  Valid (t-t_0) range: [0.00, {(t_valid.max() - t0):.2f}]")

        # Collect SHIFTED times for global fit
        t_shifted = t_valid - t0
        all_valid_t_shifted.extend(t_shifted)
        all_valid_xf.extend(xf_valid)

# Global fit: log(xf) = log(A) + 1.5 * log(t - t_0)
all_valid_t_shifted = np.array(all_valid_t_shifted)
all_valid_xf = np.array(all_valid_xf)

# Filter out t_shifted = 0 (would give -inf in log)
mask_positive = all_valid_t_shifted > 0
all_valid_t_shifted = all_valid_t_shifted[mask_positive]
all_valid_xf = all_valid_xf[mask_positive]

print(f"\n{'=' * 60}")
print(f"Total valid points for global fit: {len(all_valid_t_shifted)}")

log_t_shifted = np.log(all_valid_t_shifted)
log_xf = np.log(all_valid_xf)

# First: unconstrained fit to check actual slope
slope_free, intercept_free, r_free, _, _ = stats.linregress(log_t_shifted, log_xf)
A_free = np.exp(intercept_free)

print(f"\n--- Unconstrained fit: x = A*(t-t_0)^n ---")
print(f"  Slope n = {slope_free:.4f}  (theory: 1.5)")
print(f"  Prefactor A = {A_free:.4f}")
print(f"  R² = {r_free**2:.6f}")

# Second: constrained fit with slope = 1.5
# log(A) = mean(log_xf - 1.5 * log(t - t_0))
log_A = np.mean(log_xf - 1.5 * log_t_shifted)
A_fitted = np.exp(log_A)

# For constrained fit, compute R² in log-space (more appropriate)
log_xf_predicted = log_A + 1.5 * log_t_shifted
ss_res_log = np.sum((log_xf - log_xf_predicted)**2)
ss_tot_log = np.sum((log_xf - np.mean(log_xf))**2)
r_squared_log = 1 - ss_res_log / ss_tot_log

# Also compute std of residuals in log-space
residuals_log = log_xf - log_xf_predicted
std_residuals = np.std(residuals_log)

print(f"\n--- Constrained fit: x = A*(t-t_0)^(3/2) ---")
print(f"  Prefactor A = {A_fitted:.4f}")
print(f"  R² (log-space) = {r_squared_log:.6f}")
print(f"  Std of log-residuals = {std_residuals:.4f}")
print("=" * 60)

# =============================================================================
# Plot: Compensated x(t)/x_a vs t (with fit lines using t - t_0)
# =============================================================================
fig, ax = plt.subplots(figsize=(5, 4))

for idx, t_arr, end_idx, label in plot_config:
    ax.loglog(
        t_arr[1:end_idx], min_x[idx][1:end_idx] / x_a[idx],
        'o', markersize=marker_size,
        markerfacecolor=colors[idx], markeredgecolor='none',
        alpha=0.7, label=label
    )

# Add valid time window indicators (dashed vertical lines)
for idx in dataset_t_range:
    t_min, t_max = dataset_t_range[idx]
    ax.axvline(t_min, color=colors[idx], linestyle='--', linewidth=1.0, alpha=0.6)
    ax.axvline(t_max, color=colors[idx], linestyle='--', linewidth=1.0, alpha=0.6)

# Single reference line showing the scaling (using representative t_0)
t_ref = np.logspace(-1, 2.5, 400)
ax.loglog(t_ref, A_fitted * t_ref**1.5, 'k-', linewidth=2.5,
          label=rf'${A_fitted:.2f}\,t^{{3/2}}$')

ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$x/(c_v\beta\Gamma_0\sqrt{\mathrm{Pe}}\theta^4/(3\sqrt{\pi}))$')
ax.legend(loc='lower right', fontsize=8)

fig.tight_layout()
output_path = Path(__file__).parent / 'fit_prefactor.pdf'
fig.savefig(output_path)
print(f"\nSaved: {output_path}")

# =============================================================================
# Summary of t_0 values
# =============================================================================
print(f"\n{'=' * 60}")
print("Summary of t_0 values by dataset:")
print("=" * 60)
for idx, t_arr, end_idx, label in plot_config:
    if idx in dataset_t0:
        print(f"  Dataset {idx}: t_0 = {dataset_t0[idx]:.2f}  ({label})")

plt.close('all')
