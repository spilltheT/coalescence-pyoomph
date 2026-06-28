"""
Marangoni stress visualization for surfactant coalescence simulations.

Plots the surfactant distribution Γ(x) and Marangoni stress T(x) = -β∂Γ/∂x
at a specified time, and computes the cancellation index χ.

Physics:
    The Marangoni stress T = ∂γ/∂s ≈ -β ∂Γ/∂x becomes sign-changing when
    Γ(x) develops a local maximum. The cancellation index χ quantifies
    the "partial cancellation" of opposing stress lobes:

        χ = 1 - |I_T| / I_{|T|}  ∈ [0, 1]

    where I_T = ∫_B T dx (signed) and I_{|T|} = ∫_B |T| dx (unsigned).
    - χ ≈ 0: one-signed stress (no cancellation)
    - χ ≈ 1: two opposite-sign lobes (strong cancellation)

Usage:
    python plot_marangoni.py <folder> --time 10.0 --beta 0.1
    python plot_marangoni.py coalescence --time 50.0 --beta 0.1 --bridge-width 0.5

Output:
    Saves plots to {folder}_plots/:
    - marangoni_t{time}.pdf: Γ(x) and Marangoni stress T(x)
    - velocity_t{time}.pdf: Surface velocity u_s and its components
"""
import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from postprocess_functions import plt_settings, style_axis


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot Marangoni stress and compute cancellation index χ",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("folder", type=str,
                        help="Simulation output folder (e.g., coalescence)")
    parser.add_argument("--time", type=float, required=True,
                        help="Time snapshot to plot")
    parser.add_argument("--beta", type=float, required=True,
                        help="Surfactant strength (β) for T = -β ∂Γ/∂x")
    parser.add_argument("--bridge-width", type=float, default=1.0,
                        help="Half-width of bridge region for χ calculation (|x| < width)")
    parser.add_argument("--output", type=str, default=None,
                        help="Output filename (default: {folder}_plots/marangoni_t{time}.pdf)")
    parser.add_argument("--png", action="store_true",
                        help="Save as PNG instead of PDF")
    return parser.parse_args()


def find_closest_file(folder, target_time):
    """Find the data file with time closest to target_time."""
    domain_dir = os.path.join(folder, "domain")
    if not os.path.exists(domain_dir):
        print(f"Error: {domain_dir} not found")
        sys.exit(1)

    files = sorted([f for f in os.listdir(domain_dir) if f.endswith('.txt')])
    if not files:
        print(f"Error: No .txt files in {domain_dir}")
        sys.exit(1)

    best_file = None
    best_time = None
    min_diff = float('inf')

    for f in files:
        filepath = os.path.join(domain_dir, f)
        with open(filepath) as file:
            header = file.readline()
            time = float(header.split('@time=')[-1])
            diff = abs(time - target_time)
            if diff < min_diff:
                min_diff = diff
                best_file = filepath
                best_time = time

    return best_file, best_time


def load_data(filepath):
    """Load x, h, p, Gamma from data file."""
    with open(filepath) as f:
        f.readline()  # skip header
        data = np.loadtxt(f)

    x = data[:, 0]
    h = data[:, 1]
    p = data[:, 2]
    Gamma = data[:, 3]

    # Sort by x
    sort_idx = np.argsort(x)
    return x[sort_idx], h[sort_idx], p[sort_idx], Gamma[sort_idx]


def compute_chi(x, T, bridge_width):
    """Compute the cancellation index χ over the bridge region."""
    mask = np.abs(x) < bridge_width
    x_bridge = x[mask]
    T_bridge = T[mask]

    if len(x_bridge) < 2:
        return np.nan, 0.0, 0.0

    I_T = np.trapezoid(T_bridge, x_bridge)          # signed integral
    I_abs_T = np.trapezoid(np.abs(T_bridge), x_bridge)  # unsigned integral

    if I_abs_T == 0:
        return np.nan, I_T, I_abs_T

    chi = 1 - np.abs(I_T) / I_abs_T
    return chi, I_T, I_abs_T


def main():
    args = parse_args()

    # Find closest file to requested time
    filepath, actual_time = find_closest_file(args.folder, args.time)
    print(f"Loading data from: {filepath}")
    print(f"Requested time: {args.time}, Actual time: {actual_time:.4f}")

    # Load data
    x, h, p, Gamma = load_data(filepath)

    # Compute Marangoni stress: T = -β ∂Γ/∂x
    dGamma_dx = np.gradient(Gamma, x)
    T = -args.beta * dGamma_dx

    # Compute surface velocity components
    dp_dx = np.gradient(p, x)
    u_pressure = -h**2 / 2 * dp_dx           # Pressure-driven component
    u_marangoni = -args.beta * h * dGamma_dx # Marangoni component
    u_s = u_pressure + u_marangoni           # Total surface velocity

    # Compute cancellation index
    chi, I_T, I_abs_T = compute_chi(x, T, args.bridge_width)

    print(f"\nCancellation Analysis (|x| < {args.bridge_width}):")
    print(f"  I_T (signed)   = {I_T:.6f}")
    print(f"  I_|T| (unsigned) = {I_abs_T:.6f}")
    print(f"  χ = 1 - |I_T|/I_|T| = {chi:.4f}")
    if chi > 0.7:
        print("  → Strong cancellation (opposing lobes)")
    elif chi < 0.3:
        print("  → Weak cancellation (one-signed stress)")
    else:
        print("  → Moderate cancellation")

    # Create figure with 2 subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    # Plot 1: Γ(x)
    ax1.plot(x, Gamma, 'b-', linewidth=2.5)
    style_axis(ax1, ylabel=r'Surfactant $\Gamma$',
               title=rf'Surfactant and Marangoni Stress at $t = {actual_time:.2f}$')

    # Plot 2: T(x) = -β ∂Γ/∂x
    ax2.plot(x, T, 'r-', linewidth=2.5)
    ax2.axhline(y=0, color='k', linestyle='-', alpha=0.3, linewidth=1)

    # Restrict to bridge region
    ax2.set_xlim(-args.bridge_width, args.bridge_width)

    # Fill positive and negative regions differently
    mask_bridge = np.abs(x) < args.bridge_width
    x_b = x[mask_bridge]
    T_b = T[mask_bridge]
    ax2.fill_between(x_b, T_b, 0, where=(T_b > 0), alpha=0.3, color='green', label=r'$T > 0$')
    ax2.fill_between(x_b, T_b, 0, where=(T_b < 0), alpha=0.3, color='purple', label=r'$T < 0$')

    style_axis(ax2, xlabel=r'Position $x$',
               ylabel=rf'Marangoni stress $T = -\beta \partial\Gamma/\partial x$')
    ax2.legend(fontsize=plt_settings['LegendFont'], frameon=False, loc='upper right')

    # Add χ annotation
    ax2.text(0.02, 0.95, rf'$\chi = {chi:.3f}$', transform=ax2.transAxes,
             fontsize=plt_settings['LabelFont'], verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()

    # Save figure
    plot_dir = f"{args.folder}_plots"
    os.makedirs(plot_dir, exist_ok=True)
    ext = "png" if args.png else "pdf"

    if args.output:
        output_path = args.output
    else:
        output_path = os.path.join(plot_dir, f"marangoni_t{actual_time:.1f}.{ext}")

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {output_path}")
    plt.close()

    # Create velocity figure
    fig2, ax_vel = plt.subplots(figsize=(10, 5))

    # Plot all three curves
    ax_vel.plot(x, u_s, '-', color='#1f77b4', linewidth=2.5, label=r'$u_s$')
    ax_vel.plot(x, u_pressure, '-', color='#ff7f0e', linewidth=2.5,
                label=r'$-h^2 \partial_x p/2$')
    ax_vel.plot(x, u_marangoni, '-', color='#9467bd', linewidth=2.5,
                label=r'$-\beta h \partial_x \Gamma$')

    ax_vel.axhline(y=0, color='gray', linestyle='-', alpha=0.5, linewidth=1)
    ax_vel.set_xlim(-args.bridge_width, args.bridge_width)

    style_axis(ax_vel, xlabel=r'Position $x$',
               ylabel=r'Surface velocity $u_s$',
               title=rf'Surface Velocity Components at $t = {actual_time:.2f}$')
    ax_vel.legend(fontsize=plt_settings['LegendFont'], frameon=False, loc='best')

    plt.tight_layout()

    # Save velocity figure
    vel_output_path = os.path.join(plot_dir, f"velocity_t{actual_time:.1f}.{ext}")
    plt.savefig(vel_output_path, dpi=300, bbox_inches='tight')
    print(f"Velocity plot saved to: {vel_output_path}")
    plt.close()


if __name__ == "__main__":
    main()
