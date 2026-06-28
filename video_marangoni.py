#!/usr/bin/env python
"""
Generate video frames for Marangoni stress visualization.

Creates frames showing:
  - Top panel: Height h(x) and surfactant Γ(x) with dual y-axes
  - Bottom panel: Normalized Marangoni stress T/T_max with fill regions
    (T = -β ∂Γ/∂x, normalized so upper limit is always 1)

Uses two-pass approach: first scan for global T_max, then generate frames.

Usage:
    python video_marangoni.py coalescence --beta 0.1
    python video_marangoni.py sweep_Pe_beta0.1/Pe_100 --beta 0.1 --bridge-width 0.5

Output:
    Frames saved to {folder}_frames-marangoni/frame_XXXXX.png
"""
import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from postprocess_functions import plt_settings, style_axis


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate video frames for Marangoni stress visualization",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("folder", type=str,
                        help="Simulation output folder (e.g., coalescence)")
    parser.add_argument("--beta", type=float, required=True,
                        help="Surfactant strength (β) for T = -β ∂Γ/∂x")
    parser.add_argument("--bridge-width", type=float, default=1.0,
                        help="Half-width of x-axis region (|x| < width)")
    return parser.parse_args()


def load_data(filepath):
    """Load x, h, Gamma from data file and return with time."""
    with open(filepath) as f:
        header = f.readline()
        time = float(header.split('@time=')[-1])
        data = np.loadtxt(f)

    x = data[:, 0]
    h = data[:, 1]
    Gamma = data[:, 3]

    # Sort by x
    sort_idx = np.argsort(x)
    return time, x[sort_idx], h[sort_idx], Gamma[sort_idx]


def compute_limits(folder, beta, bridge_width):
    """
    First pass: scan all files to determine global min/max values.
    Returns dict with limits for h, Gamma, and T_max for normalization.
    """
    domain_dir = os.path.join(folder, "domain")
    files = sorted([f for f in os.listdir(domain_dir) if f.endswith('.txt')])

    h_min, h_max = float('inf'), float('-inf')
    Gamma_min, Gamma_max = float('inf'), float('-inf')
    T_max = float('-inf')

    # Skip first few files for T limits (initial step function causes spike)
    skip_initial = 5

    print(f"Scanning {len(files)} files for axis limits...")

    for i, f in enumerate(files):
        filepath = os.path.join(domain_dir, f)
        _, x, h, Gamma = load_data(filepath)

        # Restrict to bridge region
        mask = np.abs(x) < bridge_width
        x_b, h_b, Gamma_b = x[mask], h[mask], Gamma[mask]

        # Compute Marangoni stress
        dGamma_dx = np.gradient(Gamma_b, x_b)
        T_b = -beta * dGamma_dx

        # Update limits (skip initial files for T due to step function spike)
        h_min = min(h_min, np.min(h_b))
        h_max = max(h_max, np.max(h_b))
        Gamma_min = min(Gamma_min, np.min(Gamma_b))
        Gamma_max = max(Gamma_max, np.max(Gamma_b))
        if i >= skip_initial:
            T_max = max(T_max, np.max(T_b))

        if (i + 1) % 100 == 0:
            print(f"  Scanned {i + 1}/{len(files)} files...")

    # Add some padding
    h_pad = 0.05 * (h_max - h_min)
    Gamma_pad = 0.05 * (Gamma_max - Gamma_min)

    limits = {
        'h': (h_min - h_pad, h_max + h_pad),
        'Gamma': (Gamma_min - Gamma_pad, Gamma_max + Gamma_pad),
        'T_max': T_max,
    }

    print(f"  h range: [{limits['h'][0]:.4f}, {limits['h'][1]:.4f}]")
    print(f"  Γ range: [{limits['Gamma'][0]:.4f}, {limits['Gamma'][1]:.4f}]")
    print(f"  T_max: {limits['T_max']:.4f} (normalizing T/T_max, y-limits: -0.2 to 1.0)")

    return limits


def generate_frame(filepath, limits, beta, bridge_width, output_path):
    """Generate a single frame with dual-axis top plot and T(x) bottom plot."""
    time, x, h, Gamma = load_data(filepath)

    # Restrict to bridge region
    mask = np.abs(x) < bridge_width
    x_b, h_b, Gamma_b = x[mask], h[mask], Gamma[mask]

    # Compute Marangoni stress and normalize by T_max
    dGamma_dx = np.gradient(Gamma_b, x_b)
    T_b = -beta * dGamma_dx
    T_normalized = T_b / limits['T_max']

    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

    # ===== Top subplot: h(x) and Γ(x) with dual y-axes =====
    color_h = '#1f77b4'  # blue
    color_Gamma = '#d62728'  # red

    # Height on left axis
    line1, = ax1.plot(x_b, h_b, color=color_h, linewidth=2.5, label=r'$h(x)$')
    ax1.set_ylabel(r'Height $h$', color=color_h, fontsize=plt_settings['LabelFont'])
    ax1.tick_params(axis='y', labelcolor=color_h, labelsize=plt_settings['AxesFont'])
    ax1.set_ylim(limits['h'])

    # Surfactant on right axis
    ax1_twin = ax1.twinx()
    line2, = ax1_twin.plot(x_b, Gamma_b, color=color_Gamma, linewidth=2.5, label=r'$\Gamma(x)$')
    ax1_twin.set_ylabel(r'Surfactant $\Gamma$', color=color_Gamma, fontsize=plt_settings['LabelFont'])
    ax1_twin.tick_params(axis='y', labelcolor=color_Gamma, labelsize=plt_settings['AxesFont'])
    ax1_twin.set_ylim(limits['Gamma'])

    # X-axis and title
    ax1.set_xlim(-bridge_width, bridge_width)
    ax1.set_xlabel(r'Position $x$', fontsize=plt_settings['LabelFont'])
    ax1.tick_params(axis='x', labelsize=plt_settings['AxesFont'])
    ax1.set_title(rf'$t = {time:.2f}$', fontsize=plt_settings['TitleFont'], pad=10)

    # Combined legend
    lines = [line1, line2]
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='upper right', fontsize=plt_settings['LegendFont'], frameon=False)

    # Styling for spines
    for spine in ax1.spines.values():
        spine.set_linewidth(2)
    for spine in ax1_twin.spines.values():
        spine.set_linewidth(2)

    # ===== Bottom subplot: Normalized Marangoni stress T/T_max =====
    ax2.plot(x_b, T_normalized, 'k-', linewidth=2.5)
    ax2.axhline(y=0, color='gray', linestyle='-', alpha=0.5, linewidth=1)

    # Fill positive and negative regions
    ax2.fill_between(x_b, T_normalized, 0, where=(T_normalized > 0), alpha=0.4, color='green', label=r'$T > 0$')
    ax2.fill_between(x_b, T_normalized, 0, where=(T_normalized < 0), alpha=0.4, color='purple', label=r'$T < 0$')

    ax2.set_xlim(-bridge_width, bridge_width)
    ax2.set_ylim(-0.2, 1.0)
    style_axis(ax2, xlabel=r'Position $x$',
               ylabel=r'Normalized stress $T/T_{\max}$')
    ax2.legend(fontsize=plt_settings['LegendFont'], frameon=False, loc='upper right')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def main():
    args = parse_args()

    # Validate folder
    domain_dir = os.path.join(args.folder, "domain")
    if not os.path.exists(domain_dir):
        print(f"Error: {domain_dir} not found")
        sys.exit(1)

    files = sorted([f for f in os.listdir(domain_dir) if f.endswith('.txt')])
    if not files:
        print(f"Error: No .txt files in {domain_dir}")
        sys.exit(1)

    print(f"Processing {len(files)} files from {domain_dir}")
    print(f"Parameters: beta = {args.beta}, bridge_width = {args.bridge_width}")
    print()

    # Create output directory
    frame_dir = f"{args.folder}_frames-marangoni"
    os.makedirs(frame_dir, exist_ok=True)

    # First pass: compute global limits
    limits = compute_limits(args.folder, args.beta, args.bridge_width)
    print()

    # Second pass: generate frames
    print(f"Generating {len(files)} frames...")
    for i, f in enumerate(files):
        filepath = os.path.join(domain_dir, f)
        output_path = os.path.join(frame_dir, f'frame_{i:05d}.png')
        generate_frame(filepath, limits, args.beta, args.bridge_width, output_path)

        if (i + 1) % 50 == 0:
            print(f"  Generated {i + 1}/{len(files)} frames...")

    print(f"\nAll {len(files)} frames saved to '{frame_dir}/'")
    print()
    print("To create video:")
    print(f"  ffmpeg -framerate 30 -i '{frame_dir}/frame_%05d.png' \\")
    print(f"    -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' \\")
    print(f"    -c:v libx264 -pix_fmt yuv420p {args.folder}_marangoni.mp4")


if __name__ == "__main__":
    main()
