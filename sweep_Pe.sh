#!/bin/bash
# Parameter sweep for Péclet number
# Runs 16 cases in batches of 4: Pe = 10, 25, 50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000
#
# Usage: ./sweep_Pe.sh --beta VALUE [--process]
#   --beta VALUE  Surfactant strength (REQUIRED)
#   --process     Skip simulations, only run post-processing

set -e  # Exit on error

# Check for required --beta argument
BETA=""
PROCESS_ONLY=false
while [[ $# -gt 0 ]]; do
    case $1 in
        --beta)
            BETA="$2"
            shift 2
            ;;
        --process)
            PROCESS_ONLY=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 --beta VALUE [--process]"
            exit 1
            ;;
    esac
done

if [[ -z "$BETA" ]]; then
    echo "Error: --beta is required"
    echo "Usage: $0 --beta VALUE [--process]"
    exit 1
fi

echo "=============================================="
echo "Parameter Sweep: Péclet Number"
echo "Surfactant strength: beta=$BETA"
if $PROCESS_ONLY; then
    echo "Mode: Post-processing only"
fi
echo "=============================================="

# Péclet number values to test (16 values)
PE_VALUES=(10 25 50 100 150 200 250 300 350 400 500 600 700 800 900 1000)
BATCH_SIZE=4

# Output base directory
BASE_DIR="sweep_Pe_beta${BETA}"
PLOT_DIR="${BASE_DIR}/plots"

if $PROCESS_ONLY; then
    # Check that simulation data exists
    if [[ ! -d "$BASE_DIR" ]]; then
        echo "Error: $BASE_DIR not found. Run simulations first without --process flag."
        exit 1
    fi
    mkdir -p "$PLOT_DIR"
else
    # Clean up old output directories
    echo ""
    echo "Cleaning up old output directories..."
    rm -rf "$BASE_DIR"
    mkdir -p "$BASE_DIR"
    mkdir -p "$PLOT_DIR"

    # Launch simulations in batches
    echo ""
    echo "Launching ${#PE_VALUES[@]} simulations in batches of $BATCH_SIZE..."
    echo "----------------------------------------------"

    for ((i=0; i<${#PE_VALUES[@]}; i+=BATCH_SIZE)); do
        PIDS=()
        BATCH_PE=()

        # Launch batch
        for ((j=i; j<i+BATCH_SIZE && j<${#PE_VALUES[@]}; j++)); do
            PE=${PE_VALUES[$j]}
            OUTDIR="${BASE_DIR}/Pe_${PE}"
            LOGFILE="${BASE_DIR}/Pe_${PE}_log.txt"

            python coalescence.py \
                --beta $BETA \
                --Pe $PE \
                --output-dir "$OUTDIR" > "$LOGFILE" 2>&1 &
            PIDS+=($!)
            BATCH_PE+=($PE)
            echo "Started Pe = $PE (PID: $!)"
        done

        # Wait for batch to complete
        echo "Waiting for batch (Pe = ${BATCH_PE[*]})..."
        for k in "${!PIDS[@]}"; do
            wait ${PIDS[$k]}
            echo "Completed Pe = ${BATCH_PE[$k]}"
        done
        echo ""
    done

    echo "All simulations completed!"
fi

echo ""
echo "Generating comparison plots..."

# Create comparison plots using Python
export BETA
python << 'PYTHON_SCRIPT'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
from postprocess_functions import plt_settings, style_axis, find_neck, find_drop_edge

# Parameters
beta = float(os.environ['BETA'])
pe_values = [10, 25, 50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000]
base_dir = f'sweep_Pe_beta{beta}'
plot_dir = f'{base_dir}/plots'

# Use log scale for colormap to handle wide range of Pe values
log_pe = np.log10(pe_values)
norm = plt.Normalize(log_pe.min(), log_pe.max())
cmap = cm.viridis

# Extract time series data for each simulation
all_data = {}
hp = 1e-4  # default precursor film thickness

for pe in pe_values:
    domain_dir = f'{base_dir}/Pe_{pe}/domain'
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

    all_data[pe] = {
        'time': np.array(time_data), 'h0': np.array(h0_data),
        'x0': np.array(x0_data), 'xe': np.array(xe_data),
        'color': cmap(norm(np.log10(pe)))
    }

print("\n" + "="*60)
print(f"PARAMETER SWEEP: PÉCLET NUMBER (beta = {beta})")
print("="*60)

# Print summary for each Pe
for pe in pe_values:
    if pe in all_data:
        data = all_data[pe]
        print(f"Pe = {pe:4d}: Final h0 = {data['h0'][-1]:.6f}, x0 = {data['x0'][-1]:.6f}, xe = {data['xe'][-1]:.6f}")

# Generate individual plots
fig1, ax1 = plt.subplots(figsize=(10, 8))
for pe in pe_values:
    if pe in all_data:
        data = all_data[pe]
        ax1.plot(data['time'], data['h0'], color=data['color'], linewidth=1.5)
style_axis(ax1, xlabel=r'Time $t$', ylabel=r'Neck height $h_0$',
           title=rf'Pe Sweep ($\beta = {beta}$): $h_0(t)$')
ax1.set_xlim(left=0)
# Add colorbar
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax1)
cbar.set_label(r'$\log_{10}(\mathrm{Pe})$', fontsize=plt_settings['LabelFont'])
cbar.ax.tick_params(labelsize=plt_settings['AxesFont'])
plt.tight_layout()
plt.savefig(f'{plot_dir}/h0_vs_time.pdf', dpi=300, bbox_inches='tight')
plt.close()

fig2, ax2 = plt.subplots(figsize=(10, 8))
for pe in pe_values:
    if pe in all_data:
        data = all_data[pe]
        ax2.plot(data['time'], data['x0'], color=data['color'], linewidth=1.5)
style_axis(ax2, xlabel=r'Time $t$', ylabel=r'Neck position $x_0$',
           title=rf'Pe Sweep ($\beta = {beta}$): $x_0(t)$')
ax2.set_xlim(left=0)
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax2)
cbar.set_label(r'$\log_{10}(\mathrm{Pe})$', fontsize=plt_settings['LabelFont'])
cbar.ax.tick_params(labelsize=plt_settings['AxesFont'])
plt.tight_layout()
plt.savefig(f'{plot_dir}/x0_vs_time.pdf', dpi=300, bbox_inches='tight')
plt.close()

fig3, ax3 = plt.subplots(figsize=(10, 8))
for pe in pe_values:
    if pe in all_data:
        data = all_data[pe]
        ax3.plot(data['time'], data['xe'], color=data['color'], linewidth=1.5)
style_axis(ax3, xlabel=r'Time $t$', ylabel=r'Drop edge $x_e$',
           title=rf'Pe Sweep ($\beta = {beta}$): $x_e(t)$')
ax3.set_xlim(left=0)
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax3)
cbar.set_label(r'$\log_{10}(\mathrm{Pe})$', fontsize=plt_settings['LabelFont'])
cbar.ax.tick_params(labelsize=plt_settings['AxesFont'])
plt.tight_layout()
plt.savefig(f'{plot_dir}/xe_vs_time.pdf', dpi=300, bbox_inches='tight')
plt.close()

# Combined 3-panel figure
fig4, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 14), sharex=True)
for pe in pe_values:
    if pe in all_data:
        data = all_data[pe]
        ax1.plot(data['time'], data['h0'], color=data['color'], linewidth=1.5)
        ax2.plot(data['time'], data['x0'], color=data['color'], linewidth=1.5)
        ax3.plot(data['time'], data['xe'], color=data['color'], linewidth=1.5)

style_axis(ax1, ylabel=r'$h_0$', title=rf'Péclet Number Sweep ($\beta = {beta}$)')
style_axis(ax2, ylabel=r'$x_0$')
style_axis(ax3, xlabel=r'Time $t$', ylabel=r'$x_e$')
ax3.set_xlim(left=0)

# Add colorbar on the right side
fig4.subplots_adjust(right=0.85)
cbar_ax = fig4.add_axes([0.88, 0.15, 0.02, 0.7])
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig4.colorbar(sm, cax=cbar_ax)
cbar.set_label(r'$\log_{10}(\mathrm{Pe})$', fontsize=plt_settings['LabelFont'])
cbar.ax.tick_params(labelsize=plt_settings['AxesFont'])

plt.savefig(f'{plot_dir}/sweep_Pe_combined.pdf', dpi=300, bbox_inches='tight')
plt.close()

print(f"\nPlots saved to '{plot_dir}/' directory")
PYTHON_SCRIPT

# Run postprocessing for each individual case
echo ""
echo "Running postprocessing for each case..."
for PE in "${PE_VALUES[@]}"; do
    OUTDIR="${BASE_DIR}/Pe_${PE}"
    if [[ -d "$OUTDIR" ]]; then
        echo "  Postprocessing $OUTDIR..."
        python postprocess.py "$OUTDIR"
    fi
done

echo ""
echo "Parameter sweep complete!"
echo "Results saved to:"
echo "  - ${PLOT_DIR}/ (comparison plots)"
for PE in "${PE_VALUES[@]}"; do
    echo "  - ${BASE_DIR}/Pe_${PE}_plots/ (individual case plots)"
done
