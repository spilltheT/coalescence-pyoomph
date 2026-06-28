"""
Postprocessing for clean coalescence and spreading simulations.

Handles both clean drop coalescence (no surfactant) and axisymmetric
droplet spreading cases. Auto-detects mode based on directory name.

Usage:
    python postprocess_clean.py <output_directory>
    python postprocess_clean.py coalescence_clean  # clean coalescence
    python postprocess_clean.py spreading_clean    # axisymmetric spreading

Output Files (saved to <directory>_plots/):

    For both modes:
        1. height_pressure_evolution.pdf
           - Top: height profiles h(x) at selected times
           - Bottom: pressure profiles p(x) at same times

        2. spacetime_height.pdf
           - Contour plot of h(x,t) or h(r,t)

        3. time_series_analysis.pdf
           - Coalescence: h_0(t), x_0(t), h_max(t), p_max(t)
           - Spreading: h_max(t), h_center(t), p_max(t)

        4. zoom_region.pdf
           - Coalescence: bridge region |x| < 0.5
           - Spreading: contact line region r < 2

    Coalescence only:
        5. neck_trajectory.pdf - phase portrait (x_0, h_0)

    Spreading only:
        5. phase_portrait.pdf - (h_max, p_max) trajectory
        6. contact_line_position.pdf - contact radius r_c(t)

Mode Detection:
    - If "spreading" appears in directory name: axisymmetric mode
    - Otherwise: Cartesian coalescence mode
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from matplotlib.gridspec import GridSpec
from postprocess_functions import plt_settings, style_axis, find_neck, find_drop_edge

# Get base folder from command line, default to "coalescence_clean"
base_folder = sys.argv[1] if len(sys.argv) > 1 else "coalescence_clean"
output_dir = f"{base_folder}/domain"
plot_dir = f"{base_folder}_plots"

# Detect if this is a spreading (axisymmetric) or coalescence case
is_spreading = "spreading" in base_folder.lower()

# Check if output directory exists before listing files
if not os.path.isdir(output_dir):
    print(f"Error: Output directory not found: {output_dir}")
    print(f"Please run the simulation first to generate output files.")
    sys.exit(1)

files = sorted([f for f in os.listdir(output_dir) if f.endswith('.txt')])

# Create folder for plots
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

print(f"Processing {len(files)} files from {output_dir}...")
print(f"Mode: {'Spreading (axisymmetric)' if is_spreading else 'Coalescence'}")

# Initialize lists for time series data
time_data = []
h0_data = []       # neck height (coalescence) or max height (spreading)
x0_data = []       # neck position (coalescence only)
h_max_data = []
h_center_data = []
p_max_data = []

# Plot 1: Height profile evolution at selected times
print("Creating height evolution plot...")
fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True, constrained_layout=True)

# Select files to plot (every 40th file, or fewer if not enough files)
step = max(1, len(files) // 10)
selected_indices = range(0, len(files), step)
colors = plt.cm.viridis(np.linspace(0, 1, len(list(selected_indices))))

time_values = []
for idx, file_idx in enumerate(range(0, len(files), step)):
    with open(os.path.join(output_dir, files[file_idx])) as f:
        header = f.readline()
        time = float(header.split('@time=')[-1])
        time_values.append(time)
        data = np.loadtxt(f)
        x = data[:, 0]
        h = data[:, 1]
        p = data[:, 2]

        # Sort by x for clean plotting
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        h = h[sort_idx]
        p = p[sort_idx]

        ax1.plot(x, h, color=colors[idx])
        ax2.plot(x, p, color=colors[idx])

title = r'Droplet Spreading Evolution' if is_spreading else r'Droplet Coalescence Evolution (Clean)'
style_axis(ax1, ylabel=r'Height $h$', title=title)

xlabel = r'Radial position $r$' if is_spreading else r'Position $x$'
style_axis(ax2, xlabel=xlabel, ylabel=r'Pressure $p$')

# Add colorbar for time
sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=time_values[0], vmax=time_values[-1]))
sm.set_array([])
cbar = fig1.colorbar(sm, ax=[ax1, ax2], location='right', shrink=0.6, pad=0.08)
cbar.set_label(r'Time $t$', fontsize=plt_settings['ColorbarFont'], labelpad=10)
cbar.ax.tick_params(labelsize=plt_settings['AxesFont'])

plt.savefig(f'{plot_dir}/height_pressure_evolution.pdf', dpi=300)
plt.close()

# Plot 2: Spacetime diagram of height
print("Creating spacetime diagram...")
times = []
height_matrix = []

# First pass: determine common x grid from first file
with open(os.path.join(output_dir, files[0])) as file:
    file.readline()
    data = np.loadtxt(file)
    x_data = data[:, 0]
    sort_idx = np.argsort(x_data)
    x_sorted = x_data[sort_idx]
    x_common = np.linspace(x_sorted.min(), x_sorted.max(), 500)

for i, f in enumerate(files[::2]):  # Use every other file for speed
    with open(os.path.join(output_dir, f)) as file:
        header = file.readline()
        time = float(header.split('@time=')[-1])
        times.append(time)
        data = np.loadtxt(file)
        x_data = data[:, 0]
        h_data = data[:, 1]
        # Sort and interpolate onto common grid
        sort_idx = np.argsort(x_data)
        h_interp = np.interp(x_common, x_data[sort_idx], h_data[sort_idx])
        height_matrix.append(h_interp)

height_matrix = np.array(height_matrix)

fig2, ax = plt.subplots(figsize=(12, 8))
im = ax.pcolormesh(x_common, times, height_matrix, shading='auto', cmap='viridis')
xlabel = r'Radial position $r$' if is_spreading else r'Position $x$'
style_axis(ax, xlabel=xlabel, ylabel=r'Time $t$', title=r'Spacetime Evolution of Film Height')
cbar = plt.colorbar(im, ax=ax)
cbar.set_label(r'Height $h$', fontsize=plt_settings['ColorbarFont'], labelpad=10)
cbar.ax.tick_params(labelsize=plt_settings['AxesFont'])
plt.savefig(f'{plot_dir}/spacetime_height.pdf', dpi=300, bbox_inches='tight')
plt.close()

# Plot 3: Time series analysis
print("Analyzing time series data...")
for f in files:
    with open(os.path.join(output_dir, f)) as file:
        header = file.readline()
        time = float(header.split('@time=')[-1])
        data = np.loadtxt(file)
        x = data[:, 0]
        h = data[:, 1]
        p = data[:, 2]

        # Sort by x
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        h = h[sort_idx]
        p = p[sort_idx]

        time_data.append(time)
        h_max_data.append(np.max(h))
        p_max_data.append(np.max(np.abs(p)))

        if is_spreading:
            # For spreading: track max height and center height
            h0_data.append(np.max(h))
            x0_data.append(0.0)  # not used for spreading
            # Interpolate to get height at r=0 (x is already sorted)
            h_center_data.append(np.interp(0.0, x, h))
        else:
            # For coalescence: use find_neck for robust neck detection
            x0, h0 = find_neck(x, h)
            x0_data.append(x0)
            h0_data.append(h0)
            center_idx = np.argmin(np.abs(x))
            h_center_data.append(h[center_idx])

fig3 = plt.figure(figsize=(14, 10))
gs = GridSpec(2, 2, figure=fig3)

if is_spreading:
    # Spreading mode: original layout
    ax1 = fig3.add_subplot(gs[0, :])
    ax1.plot(time_data, h_max_data, 'b-', linewidth=3, label='Maximum height (center)')
    style_axis(ax1, xlabel=r'Time $t$', ylabel=r'Maximum height',
               title=r'Evolution of Droplet Height During Spreading')
    ax1.legend(fontsize=plt_settings['LegendFont'], frameon=False)

    ax2 = fig3.add_subplot(gs[1, 0])
    ax2.plot(time_data, h_max_data, 'r-', linewidth=3, label='Maximum height')
    ax2.plot(time_data, h_center_data, 'g--', linewidth=3, label='Center height')
    style_axis(ax2, xlabel=r'Time $t$', ylabel=r'Height',
               title=r'Height Evolution at Different Positions')
    ax2.legend(fontsize=plt_settings['LegendFont'], frameon=False)

    ax3 = fig3.add_subplot(gs[1, 1])
    ax3.plot(time_data, p_max_data, 'm-', linewidth=3)
    style_axis(ax3, xlabel=r'Time $t$', ylabel=r'Max $|p|$',
               title=r'Maximum Pressure')
else:
    # Coalescence mode: 2x2 grid with h0, x0, h_max, p_max
    ax1 = fig3.add_subplot(gs[0, 0])
    ax1.plot(time_data, h0_data, 'b-', linewidth=3)
    style_axis(ax1, xlabel=r'Time $t$', ylabel=r'Neck height $h_0$',
               title=r'Neck Height Evolution')

    ax2 = fig3.add_subplot(gs[0, 1])
    ax2.plot(time_data, x0_data, 'g-', linewidth=3)
    style_axis(ax2, xlabel=r'Time $t$', ylabel=r'Neck position $x_0$',
               title=r'Neck Position Evolution')

    ax3 = fig3.add_subplot(gs[1, 0])
    ax3.plot(time_data, h_max_data, 'r-', linewidth=3)
    style_axis(ax3, xlabel=r'Time $t$', ylabel=r'Maximum height $h_{\max}$',
               title=r'Maximum Height Evolution')

    ax4 = fig3.add_subplot(gs[1, 1])
    ax4.plot(time_data, p_max_data, 'm-', linewidth=3)
    style_axis(ax4, xlabel=r'Time $t$', ylabel=r'Max $|p|$',
               title=r'Maximum Pressure')

plt.tight_layout()
plt.savefig(f'{plot_dir}/time_series_analysis.pdf', dpi=300, bbox_inches='tight')
plt.close()

# Plot 4: Zoom on region of interest
print("Creating zoom plot...")
fig4, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

# Plot at different stages
stages = [0, len(files)//4, len(files)//2, 3*len(files)//4, len(files)-1]
colors = ['blue', 'green', 'orange', 'red', 'purple']
labels = ['Initial', 'Early', 'Middle', 'Late', 'Final']

for idx, (stage, color, label) in enumerate(zip(stages, colors, labels)):
    with open(os.path.join(output_dir, files[stage])) as f:
        header = f.readline()
        time = float(header.split('@time=')[-1])
        data = np.loadtxt(f)
        x = data[:, 0]
        h = data[:, 1]
        p = data[:, 2]

        # Sort by x
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        h = h[sort_idx]
        p = p[sort_idx]

        if is_spreading:
            # For spreading: show full profile, focus on contact line
            mask = x < 2.0  # zoom on inner region
        else:
            # For coalescence: zoom on bridge at x=0
            mask = np.abs(x) < 0.5

        ax1.plot(x[mask], h[mask], color=color, linewidth=3, label=f'{label} ($t={time:.2f}$)')
        ax2.plot(x[mask], p[mask], color=color, linewidth=3)

if is_spreading:
    style_axis(ax1, ylabel=r'Height $h$', title=r'Contact Line Region Evolution')
else:
    style_axis(ax1, ylabel=r'Height $h$', title=r'Bridge Region Evolution During Coalescence')
ax1.legend(fontsize=plt_settings['LegendFont'], frameon=False)

xlabel = r'Radial position $r$' if is_spreading else r'Position $x$'
style_axis(ax2, xlabel=xlabel, ylabel=r'Pressure $p$')

plt.tight_layout()
plt.savefig(f'{plot_dir}/zoom_region.pdf', dpi=300, bbox_inches='tight')
plt.close()

# Plot 5: Phase portrait (height dynamics)
print("Creating phase portrait...")
fig5, ax = plt.subplots(figsize=(10, 10))

# Use color gradient for time
colors = plt.cm.plasma(np.linspace(0, 1, len(time_data)))

if is_spreading:
    # For spreading: max height vs max pressure
    ax.scatter(h_max_data, p_max_data, c=colors, alpha=0.7, s=80, edgecolors='w', linewidth=0.5, zorder=3)
    ax.plot(h_max_data, p_max_data, 'k-', alpha=0.3, linewidth=1, zorder=2)
    style_axis(ax, xlabel=r'Maximum height', ylabel=r'Maximum $|p|$',
               title=r'Phase Portrait: Height vs Pressure')
    output_name = 'phase_portrait.pdf'
else:
    # For coalescence: neck trajectory (x0 vs h0)
    ax.scatter(x0_data, h0_data, c=colors, alpha=0.7, s=80, edgecolors='w', linewidth=0.5, zorder=3)
    ax.plot(x0_data, h0_data, 'k-', alpha=0.3, linewidth=1, zorder=2)

    # Add arrows to show direction
    arrow_indices = np.linspace(0, len(x0_data)-2, 10, dtype=int)
    for idx in arrow_indices:
        ax.annotate('', xy=(x0_data[idx+1], h0_data[idx+1]),
                    xytext=(x0_data[idx], h0_data[idx]),
                    arrowprops=dict(arrowstyle='->', color='black', alpha=0.5, lw=2))

    style_axis(ax, xlabel=r'Neck position $x_0$', ylabel=r'Neck height $h_0$',
               title=r'Neck Trajectory During Coalescence')
    output_name = 'neck_trajectory.pdf'

# Add colorbar for time
sm = plt.cm.ScalarMappable(cmap=plt.cm.plasma, norm=plt.Normalize(vmin=time_data[0], vmax=time_data[-1]))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label(r'Time $t$', fontsize=plt_settings['ColorbarFont'], labelpad=10)
cbar.ax.tick_params(labelsize=plt_settings['AxesFont'])

plt.savefig(f'{plot_dir}/{output_name}', dpi=300, bbox_inches='tight')
plt.close()

# Plot 6: Contact line / bridge position over time (for spreading)
if is_spreading:
    print("Tracking contact line position...")
    contact_positions = []
    for f in files:
        with open(os.path.join(output_dir, f)) as file:
            file.readline()
            data = np.loadtxt(file)
            x = data[:, 0]
            h = data[:, 1]
            sort_idx = np.argsort(x)
            x = x[sort_idx]
            h = h[sort_idx]
            # Find contact line using shared function
            hp_default = 1e-4  # default precursor film thickness
            contact_positions.append(find_drop_edge(x, h, hp_default, side='right'))

    fig6, ax = plt.subplots(figsize=(12, 8))
    ax.plot(time_data, contact_positions, 'b-', linewidth=3)
    style_axis(ax, xlabel=r'Time $t$', ylabel=r'Contact line position $r_c$',
               title=r'Contact Line Evolution During Spreading')
    plt.savefig(f'{plot_dir}/contact_line_position.pdf', dpi=300, bbox_inches='tight')
    plt.close()

print(f"\nAll plots have been saved to the '{plot_dir}' folder!")
print(f"Total number of timesteps analyzed: {len(files)}")
print(f"Time range: {time_data[0]:.2f} to {time_data[-1]:.2f}")
if is_spreading:
    print(f"Maximum height: {max(h_max_data):.6f}")
else:
    print(f"Final neck height: {h0_data[-1]:.6f}")
    print(f"Final neck position: {x0_data[-1]:.6f}")
