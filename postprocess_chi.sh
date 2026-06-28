#!/bin/bash
# Post-process Marangoni stress cancellation index χ(t) over time
# Computes χ at each time step and generates time series plot
#
# Usage:
#   ./postprocess_chi.sh --beta VALUE [--bridge-width VALUE]
#       Process all Pe folders in sweep_Pe_beta{VALUE}/
#
#   ./postprocess_chi.sh <folder> --beta VALUE [--bridge-width VALUE]
#       Process a single simulation folder

set -e  # Exit on error

# Parse arguments
FOLDER=""
BETA=""
BRIDGE_WIDTH="1.0"

while [[ $# -gt 0 ]]; do
    case $1 in
        --beta)
            BETA="$2"
            shift 2
            ;;
        --bridge-width)
            BRIDGE_WIDTH="$2"
            shift 2
            ;;
        -*)
            echo "Unknown option: $1"
            echo "Usage: $0 [folder] --beta VALUE [--bridge-width VALUE]"
            exit 1
            ;;
        *)
            if [[ -z "$FOLDER" ]]; then
                FOLDER="$1"
            else
                echo "Error: Multiple folders specified"
                exit 1
            fi
            shift
            ;;
    esac
done

# Validate required arguments
if [[ -z "$BETA" ]]; then
    echo "Error: --beta is required"
    echo "Usage: $0 [folder] --beta VALUE [--bridge-width VALUE]"
    exit 1
fi

# Determine mode: sweep or single folder
if [[ -z "$FOLDER" ]]; then
    # Sweep mode: process all Pe folders
    SWEEP_DIR="sweep_Pe_beta${BETA}"
    if [[ ! -d "$SWEEP_DIR" ]]; then
        echo "Error: $SWEEP_DIR not found"
        exit 1
    fi

    echo "=============================================="
    echo "Marangoni Cancellation Analysis: χ(t) vs Pe"
    echo "Sweep directory: $SWEEP_DIR"
    echo "beta = $BETA, bridge_width = $BRIDGE_WIDTH"
    echo "=============================================="

    PLOT_DIR="${SWEEP_DIR}/plots"
    mkdir -p "$PLOT_DIR"

    # Process all Pe folders using Python
    export BETA BRIDGE_WIDTH SWEEP_DIR PLOT_DIR
    python << 'PYTHON_SCRIPT'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
from postprocess_functions import plt_settings, style_axis

# Parameters from environment
beta = float(os.environ['BETA'])
bridge_width = float(os.environ['BRIDGE_WIDTH'])
sweep_dir = os.environ['SWEEP_DIR']
plot_dir = os.environ['PLOT_DIR']

pe_values = [10, 25, 50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000]

# Use log scale for colormap
log_pe = np.log10(pe_values)
norm = plt.Normalize(log_pe.min(), log_pe.max())
cmap = cm.viridis

all_data = {}

for pe in pe_values:
    domain_dir = os.path.join(sweep_dir, f'Pe_{pe}', 'domain')
    if not os.path.exists(domain_dir):
        print(f"Warning: {domain_dir} not found, skipping")
        continue

    files = sorted([f for f in os.listdir(domain_dir) if f.endswith('.txt')])
    print(f"Processing Pe = {pe} ({len(files)} time steps)...")

    time_data, chi_data, I_T_data, I_abs_T_data = [], [], [], []

    for f in files:
        filepath = os.path.join(domain_dir, f)
        with open(filepath) as file:
            header = file.readline()
            time = float(header.split('@time=')[-1])
            data = np.loadtxt(file)

        x, Gamma = data[:, 0], data[:, 3]
        sort_idx = np.argsort(x)
        x, Gamma = x[sort_idx], Gamma[sort_idx]

        # Compute Marangoni stress and chi
        dGamma_dx = np.gradient(Gamma, x)
        T = -beta * dGamma_dx

        mask = np.abs(x) < bridge_width
        x_bridge, T_bridge = x[mask], T[mask]

        if len(x_bridge) >= 2:
            I_T = np.trapezoid(T_bridge, x_bridge)
            I_abs_T = np.trapezoid(np.abs(T_bridge), x_bridge)
            chi = 1 - np.abs(I_T) / I_abs_T if I_abs_T > 0 else np.nan
        else:
            I_T, I_abs_T, chi = np.nan, np.nan, np.nan

        time_data.append(time)
        chi_data.append(chi)
        I_T_data.append(I_T)
        I_abs_T_data.append(I_abs_T)

    all_data[pe] = {
        'time': np.array(time_data),
        'chi': np.array(chi_data),
        'I_T': np.array(I_T_data),
        'I_abs_T': np.array(I_abs_T_data),
        'color': cmap(norm(np.log10(pe)))
    }

    # Save individual data
    output_file = os.path.join(plot_dir, f"chi_Pe_{pe}.txt")
    header = f"# chi(t) for Pe = {pe}, beta = {beta}, bridge_width = {bridge_width}\n"
    header += "# Columns: time, chi, I_T (signed), I_|T| (unsigned)"
    np.savetxt(output_file, np.column_stack([time_data, chi_data, I_T_data, I_abs_T_data]),
               header=header, fmt='%.6e')

# Print summary
print("\n" + "="*60)
print(f"CANCELLATION INDEX SUMMARY (beta = {beta})")
print("="*60)
for pe in pe_values:
    if pe in all_data:
        d = all_data[pe]
        print(f"Pe = {pe:4d}: Final χ = {d['chi'][-1]:.4f}, Max χ = {np.nanmax(d['chi']):.4f}")

# Generate comparison plot
fig, ax = plt.subplots(figsize=(10, 8))
for pe in pe_values:
    if pe in all_data:
        d = all_data[pe]
        ax.plot(d['time'], d['chi'], color=d['color'], linewidth=1.5)

ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax.set_ylim(0, 1)
ax.set_xlim(left=0)
style_axis(ax, xlabel=r'Time $t$', ylabel=r'Cancellation index $\chi$',
           title=rf'Marangoni Cancellation vs Pe ($\beta = {beta}$)')

# Add colorbar
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label(r'$\log_{10}(\mathrm{Pe})$', fontsize=plt_settings['LabelFont'])
cbar.ax.tick_params(labelsize=plt_settings['AxesFont'])

plt.tight_layout()
plt.savefig(os.path.join(plot_dir, "chi_vs_time.pdf"), dpi=300, bbox_inches='tight')
plt.close()
print(f"\nComparison plot saved to: {plot_dir}/chi_vs_time.pdf")

# Generate combined plot with I_T comparison
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12), sharex=True)

for pe in pe_values:
    if pe in all_data:
        d = all_data[pe]
        ax1.plot(d['time'], d['chi'], color=d['color'], linewidth=1.5)
        ax2.plot(d['time'], d['I_T'], color=d['color'], linewidth=1.5)

ax1.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax1.set_ylim(0, 1)
style_axis(ax1, ylabel=r'$\chi = |I_T|/I_{|T|}$',
           title=rf'Marangoni Cancellation Analysis ($\beta = {beta}$)')

ax2.axhline(y=0, color='k', linestyle='-', alpha=0.3, linewidth=1)
style_axis(ax2, xlabel=r'Time $t$', ylabel=r'$I_T$ (signed integral)')
ax2.set_xlim(left=0)

# Add colorbar
fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.7])
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label(r'$\log_{10}(\mathrm{Pe})$', fontsize=plt_settings['LabelFont'])
cbar.ax.tick_params(labelsize=plt_settings['AxesFont'])

plt.savefig(os.path.join(plot_dir, "chi_analysis.pdf"), dpi=300, bbox_inches='tight')
plt.close()
print(f"Combined analysis saved to: {plot_dir}/chi_analysis.pdf")
PYTHON_SCRIPT

    echo ""
    echo "Post-processing complete!"
    echo "Results saved to ${PLOT_DIR}/"

else
    # Single folder mode
    if [[ ! -d "$FOLDER/domain" ]]; then
        echo "Error: $FOLDER/domain not found"
        exit 1
    fi

    echo "=============================================="
    echo "Marangoni Stress Cancellation Analysis: χ(t)"
    echo "Folder: $FOLDER"
    echo "beta = $BETA, bridge_width = $BRIDGE_WIDTH"
    echo "=============================================="

    PLOT_DIR="${FOLDER}_plots"
    mkdir -p "$PLOT_DIR"

    # Compute χ(t) using embedded Python
    export FOLDER BETA BRIDGE_WIDTH PLOT_DIR
    python << 'PYTHON_SCRIPT'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from postprocess_functions import plt_settings, style_axis

# Parameters from environment
folder = os.environ['FOLDER']
beta = float(os.environ['BETA'])
bridge_width = float(os.environ['BRIDGE_WIDTH'])
plot_dir = os.environ['PLOT_DIR']

domain_dir = os.path.join(folder, "domain")
files = sorted([f for f in os.listdir(domain_dir) if f.endswith('.txt')])

print(f"Processing {len(files)} time steps...")

time_data, chi_data, I_T_data, I_abs_T_data = [], [], [], []

for i, f in enumerate(files):
    filepath = os.path.join(domain_dir, f)
    with open(filepath) as file:
        header = file.readline()
        time = float(header.split('@time=')[-1])
        data = np.loadtxt(file)

    x, Gamma = data[:, 0], data[:, 3]
    sort_idx = np.argsort(x)
    x, Gamma = x[sort_idx], Gamma[sort_idx]

    dGamma_dx = np.gradient(Gamma, x)
    T = -beta * dGamma_dx

    mask = np.abs(x) < bridge_width
    x_bridge, T_bridge = x[mask], T[mask]

    if len(x_bridge) >= 2:
        I_T = np.trapezoid(T_bridge, x_bridge)
        I_abs_T = np.trapezoid(np.abs(T_bridge), x_bridge)
        chi = 1 - np.abs(I_T) / I_abs_T if I_abs_T > 0 else np.nan
    else:
        I_T, I_abs_T, chi = np.nan, np.nan, np.nan

    time_data.append(time)
    chi_data.append(chi)
    I_T_data.append(I_T)
    I_abs_T_data.append(I_abs_T)

    if (i + 1) % 100 == 0:
        print(f"  Processed {i + 1}/{len(files)} time steps...")

time_data = np.array(time_data)
chi_data = np.array(chi_data)
I_T_data = np.array(I_T_data)
I_abs_T_data = np.array(I_abs_T_data)

# Save data
output_file = os.path.join(plot_dir, "chi_vs_time.txt")
header = f"# chi(t), beta = {beta}, bridge_width = {bridge_width}\n"
header += "# Columns: time, chi, I_T (signed), I_|T| (unsigned)"
np.savetxt(output_file, np.column_stack([time_data, chi_data, I_T_data, I_abs_T_data]),
           header=header, fmt='%.6e')
print(f"\nData saved to: {output_file}")

print(f"\nSummary:")
print(f"  Initial χ(t=0) = {chi_data[0]:.4f}")
print(f"  Final χ(t={time_data[-1]:.1f}) = {chi_data[-1]:.4f}")
print(f"  Min χ = {np.nanmin(chi_data):.4f} at t = {time_data[np.nanargmin(chi_data)]:.2f}")
print(f"  Max χ = {np.nanmax(chi_data):.4f} at t = {time_data[np.nanargmax(chi_data)]:.2f}")

# Generate plots
fig, ax = plt.subplots(figsize=(10, 8))
ax.plot(time_data, chi_data, 'b-', linewidth=2.5)
ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax.set_ylim(0, 1)
ax.set_xlim(left=0)
style_axis(ax, xlabel=r'Time $t$', ylabel=r'Cancellation index $\chi$',
           title=rf'Marangoni Stress Cancellation ($\beta = {beta}$)')
ax.axhspan(0, 0.3, alpha=0.1, color='red', label='Weak cancellation')
ax.axhspan(0.7, 1.0, alpha=0.1, color='green', label='Strong cancellation')
ax.legend(fontsize=plt_settings['LegendFont'], frameon=False, loc='upper right')
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, "chi_vs_time.pdf"), dpi=300, bbox_inches='tight')
plt.close()
print(f"Plot saved to: {plot_dir}/chi_vs_time.pdf")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)
ax1.plot(time_data, chi_data, 'b-', linewidth=2.5)
ax1.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax1.set_ylim(0, 1)
style_axis(ax1, ylabel=r'$\chi = |I_T|/I_{|T|}$',
           title=rf'Marangoni Cancellation Analysis ($\beta = {beta}$)')

ax2.plot(time_data, I_T_data, 'g-', linewidth=2, label=r'$I_T$ (signed)')
ax2.plot(time_data, I_abs_T_data, 'r-', linewidth=2, label=r'$I_{|T|}$ (unsigned)')
ax2.axhline(y=0, color='k', linestyle='-', alpha=0.3, linewidth=1)
style_axis(ax2, xlabel=r'Time $t$', ylabel=r'Integrated Marangoni stress')
ax2.legend(fontsize=plt_settings['LegendFont'], frameon=False)
ax2.set_xlim(left=0)
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, "chi_analysis.pdf"), dpi=300, bbox_inches='tight')
plt.close()
print(f"Combined plot saved to: {plot_dir}/chi_analysis.pdf")
PYTHON_SCRIPT

    echo ""
    echo "Post-processing complete!"
    echo "Results saved to:"
    echo "  - ${PLOT_DIR}/chi_vs_time.txt (data)"
    echo "  - ${PLOT_DIR}/chi_vs_time.pdf (plot)"
    echo "  - ${PLOT_DIR}/chi_analysis.pdf (combined analysis)"
fi
