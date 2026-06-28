#!/bin/bash
# Sensitivity analysis for domain size Lx
# Runs 4 cases in parallel: 6, 8, 10, 12
#
# Usage: ./sensitivity_Lx.sh [--Pe VALUE] [--beta VALUE]
#   --Pe VALUE    Péclet number (default: 1.0)
#   --beta VALUE  Surfactant strength (default: 0.1)

set -e  # Exit on error

# Default physical parameters (can be overridden via command line)
BETA=0.1
PE=1.0

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --Pe)
            PE="$2"
            shift 2
            ;;
        --beta)
            BETA="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--Pe VALUE] [--beta VALUE]"
            exit 1
            ;;
    esac
done

echo "=============================================="
echo "Sensitivity Analysis: Domain Size"
echo "Physical parameters: Pe=$PE, beta=$BETA"
echo "=============================================="

# Fixed parameters for this study
GAMMA0=0.8
THETA=20
HP=1e-4

# Domain size values to test (4 values)
LX_VALUES=("6" "8" "10" "12")

# Clean up old output directories
for LX in "${LX_VALUES[@]}"; do
    rm -rf "sensitivity_Lx_${LX}"
done

# Launch all simulations in parallel
echo ""
echo "Launching ${#LX_VALUES[@]} simulations in parallel..."
echo "----------------------------------------------"

PIDS=()
for LX in "${LX_VALUES[@]}"; do
    OUTDIR="sensitivity_Lx_${LX}"
    python coalescence.py \
        --beta $BETA \
        --Pe $PE \
        --Gamma0 $GAMMA0 \
        --theta $THETA \
        --hp $HP \
        --Lx $LX \
        --N 5000 \
        --output-dir "$OUTDIR" > "${OUTDIR}_log.txt" 2>&1 &
    PIDS+=($!)
    echo "Started Lx = $LX (PID: $!)"
done

# Wait for all simulations to complete
echo ""
echo "Waiting for all simulations to complete..."
for i in "${!PIDS[@]}"; do
    wait ${PIDS[$i]}
    echo "Completed Lx = ${LX_VALUES[$i]}"
done

echo ""
echo "All simulations completed!"

# Run convergence analysis
echo ""
echo "Running convergence analysis..."
python check_convergence.py --Lx --Pe $PE --beta $BETA

echo ""
echo "Generating comparison plots..."

# Create comparison plot using Python
python << 'PYTHON_SCRIPT'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from postprocess_functions import plt_settings, style_axis, find_neck, find_drop_edge

# Domain size values
Lx_values = ['6', '8', '10', '12']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
labels = [r'$L_x = 6$', r'$L_x = 8$', r'$L_x = 10$', r'$L_x = 12$']

plot_dir = 'sensitivity_Lx_plots'
os.makedirs(plot_dir, exist_ok=True)

hp = 1e-4  # Precursor film thickness

# Extract time series data for each simulation
all_data = {}

for Lx, color, label in zip(Lx_values, colors, labels):
    domain_dir = f'sensitivity_Lx_{Lx}/domain'
    if not os.path.exists(domain_dir):
        print(f"Warning: {domain_dir} not found, skipping")
        continue

    files = sorted([f for f in os.listdir(domain_dir) if f.endswith('.txt')])
    time_data, h0_data, x0_data, xe_data = [], [], [], []

    for f in files:
        with open(os.path.join(domain_dir, f)) as file:
            header = file.readline()
            time = float(header.split('@time=')[-1])
            data = np.loadtxt(file)
            x, h = data[:, 0], data[:, 1]

            x0, h0 = find_neck(x, h)
            xe = find_drop_edge(x, h, hp, side='right')

            time_data.append(time)
            h0_data.append(h0)
            x0_data.append(x0)
            xe_data.append(xe)

    all_data[Lx] = {
        'time': np.array(time_data), 'h0': np.array(h0_data),
        'x0': np.array(x0_data), 'xe': np.array(xe_data),
        'color': color, 'label': label
    }

# Print sensitivity metrics
print("\n" + "="*60)
print("SENSITIVITY ANALYSIS: DOMAIN SIZE")
print("="*60)

if '6' in all_data:
    baseline = all_data['6']
    print(f"\nBaseline: Lx = 6")
    print(f"  Final h0 = {baseline['h0'][-1]:.6f}")
    print(f"  Final x0 = {baseline['x0'][-1]:.6f}")
    print(f"  Final xe = {baseline['xe'][-1]:.6f}")

    for Lx in ['8', '10', '12']:
        if Lx in all_data:
            data = all_data[Lx]
            t_common = np.linspace(0, min(baseline['time'][-1], data['time'][-1]), 100)
            h0_base = np.interp(t_common, baseline['time'], baseline['h0'])
            h0_test = np.interp(t_common, data['time'], data['h0'])
            x0_base = np.interp(t_common, baseline['time'], baseline['x0'])
            x0_test = np.interp(t_common, data['time'], data['x0'])
            xe_base = np.interp(t_common, baseline['time'], baseline['xe'])
            xe_test = np.interp(t_common, data['time'], data['xe'])

            h0_diff = np.max(np.abs(h0_test - h0_base) / (np.abs(h0_base) + 1e-10)) * 100
            x0_diff = np.max(np.abs(x0_test - x0_base) / (np.abs(x0_base) + 1e-10)) * 100
            xe_diff = np.max(np.abs(xe_test - xe_base) / (np.abs(xe_base) + 1e-10)) * 100

            print(f"\nLx = {Lx}:")
            print(f"  Final h0 = {data['h0'][-1]:.6f}, Final x0 = {data['x0'][-1]:.6f}, Final xe = {data['xe'][-1]:.6f}")
            print(f"  Max rel. diff: h0={h0_diff:.2f}%, x0={x0_diff:.2f}%, xe={xe_diff:.2f}%")

# Generate plots
fig1, ax1 = plt.subplots(figsize=(10, 8))
for data in all_data.values():
    ax1.plot(data['time'], data['h0'], color=data['color'], label=data['label'], linewidth=2.5)
style_axis(ax1, xlabel=r'Time $t$', ylabel=r'Neck height $h_0$', title=r'Sensitivity to Domain Size: $h_0(t)$')
ax1.legend(fontsize=plt_settings['LegendFont'], frameon=False)
ax1.set_xlim(left=0)
plt.tight_layout()
plt.savefig(f'{plot_dir}/h0_vs_time.pdf', dpi=300, bbox_inches='tight')
plt.close()

fig2, ax2 = plt.subplots(figsize=(10, 8))
for data in all_data.values():
    ax2.plot(data['time'], data['x0'], color=data['color'], label=data['label'], linewidth=2.5)
style_axis(ax2, xlabel=r'Time $t$', ylabel=r'Neck position $x_0$', title=r'Sensitivity to Domain Size: $x_0(t)$')
ax2.legend(fontsize=plt_settings['LegendFont'], frameon=False)
ax2.set_xlim(left=0)
plt.tight_layout()
plt.savefig(f'{plot_dir}/x0_vs_time.pdf', dpi=300, bbox_inches='tight')
plt.close()

fig3, ax3 = plt.subplots(figsize=(10, 8))
for data in all_data.values():
    ax3.plot(data['time'], data['xe'], color=data['color'], label=data['label'], linewidth=2.5)
style_axis(ax3, xlabel=r'Time $t$', ylabel=r'Drop edge $x_e$', title=r'Sensitivity to Domain Size: $x_e(t)$')
ax3.legend(fontsize=plt_settings['LegendFont'], frameon=False)
ax3.set_xlim(left=0)
plt.tight_layout()
plt.savefig(f'{plot_dir}/xe_vs_time.pdf', dpi=300, bbox_inches='tight')
plt.close()

fig4, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 14), sharex=True)
for data in all_data.values():
    ax1.plot(data['time'], data['h0'], color=data['color'], label=data['label'], linewidth=2.5)
    ax2.plot(data['time'], data['x0'], color=data['color'], linewidth=2.5)
    ax3.plot(data['time'], data['xe'], color=data['color'], linewidth=2.5)
style_axis(ax1, ylabel=r'$h_0$', title=r'Sensitivity to Domain Size $L_x$')
ax1.legend(fontsize=plt_settings['LegendFont'], frameon=False, loc='upper left')
style_axis(ax2, ylabel=r'$x_0$')
style_axis(ax3, xlabel=r'Time $t$', ylabel=r'$x_e$')
ax3.set_xlim(left=0)
plt.tight_layout()
plt.savefig(f'{plot_dir}/sensitivity_Lx_combined.pdf', dpi=300, bbox_inches='tight')
plt.close()

print(f"\nPlots saved to '{plot_dir}/' directory")
PYTHON_SCRIPT

# Run postprocessing for each individual case
echo ""
echo "Running postprocessing for each case..."
for LX in "${LX_VALUES[@]}"; do
    OUTDIR="sensitivity_Lx_${LX}"
    echo "  Postprocessing $OUTDIR..."
    python postprocess.py "$OUTDIR"
done

echo ""
echo "Sensitivity analysis complete!"
echo "Results saved to:"
echo "  - sensitivity_Lx_plots/ (comparison plots)"
for LX in "${LX_VALUES[@]}"; do
    echo "  - sensitivity_Lx_${LX}_plots/ (individual case plots)"
done
