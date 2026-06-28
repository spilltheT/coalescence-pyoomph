"""
Postprocessing for surfactant coalescence simulations.

Generates 5 publication-quality PDF figures from simulation output,
analyzing the evolution of film height, surfactant concentration,
and coalescence dynamics.

Usage:
    python postprocess.py <output_directory>
    python postprocess.py coalescence  # default

Output Files (saved to <directory>_plots/):
    1. height_surfactant_evolution.pdf
       - Top panel: height profiles h(x) at selected times
       - Bottom panel: surfactant profiles Gamma(x) at same times
       - Time indicated by colorbar (viridis colormap)

    2. spacetime_height.pdf
       - Contour plot of h(x,t) showing coalescence progression
       - x-axis: spatial position, y-axis: time
       - Useful for visualizing bridge growth and neck motion

    3. time_series_analysis.pdf
       - 2x2 grid of key quantities vs time:
         * Neck height h_0(t): minimum film thickness
         * Neck position x_0(t): location of minimum
         * Maximum height h_max(t)
         * Maximum surfactant |Gamma|_max(t)

    4. bridge_region_zoom.pdf
       - Detailed view of bridge region (|x| < 0.5) at 5 stages
       - Shows h, p, Gamma profiles to analyze Marangoni effects
       - Stages: Initial, Early, Middle, Late, Final

    5. neck_trajectory.pdf
       - Phase portrait: x_0 vs h_0 parametrized by time
       - Arrows indicate direction of evolution
       - Colorbar shows time progression

Key Metrics Extracted:
    - h_0(t): Neck height (from find_neck using peak detection)
    - x_0(t): Horizontal position of neck minimum
    - Gamma_max(t): Maximum surfactant concentration
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from matplotlib.gridspec import GridSpec
from postprocess_functions import plt_settings, style_axis, find_neck

# Get base folder from command line, default to "coalescence"
base_folder = sys.argv[1] if len(sys.argv) > 1 else "coalescence"
output_dir = f"{base_folder}/domain"
plot_dir = f"{base_folder}_plots"

files = sorted([f for f in os.listdir(output_dir) if f.endswith('.txt')])

# Create folder for plots
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

# Initialize lists for time series data
time_data = []
h0_data = []       # neck height
x0_data = []       # neck position
h_max_data = []
gamma_max_data = []

# Read first and last files to understand the evolution
with open(os.path.join(output_dir, files[0])) as f:
    header = f.readline()
    data_first = np.loadtxt(f)
    
with open(os.path.join(output_dir, files[-1])) as f:
    header = f.readline()
    data_last = np.loadtxt(f)

# Plot 1: Height profile evolution at selected times
fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True, constrained_layout=True)

# Select files to plot (every 40th file)
selected_indices = range(0, len(files), 40)
colors = plt.cm.viridis(np.linspace(0, 1, len(selected_indices)))

time_values = []
for idx, file_idx in enumerate(selected_indices):
    with open(os.path.join(output_dir, files[file_idx])) as f:
        header = f.readline()
        time = float(header.split('@time=')[-1])
        time_values.append(time)
        data = np.loadtxt(f)
        x = data[:, 0]
        h = data[:, 1]
        gamma = data[:, 3]

        ax1.plot(x, h, color=colors[idx])
        ax2.plot(x, gamma, color=colors[idx])

style_axis(ax1, ylabel=r'Height $h$', title=r'Evolution of Droplet Coalescence with Surfactants')
style_axis(ax2, xlabel=r'Position $x$', ylabel=r'Surfactant concentration $\Gamma$')

# Add colorbar for time
sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=time_values[0], vmax=time_values[-1]))
sm.set_array([])
cbar = fig1.colorbar(sm, ax=[ax1, ax2], location='right', shrink=0.6, pad=0.08)
cbar.set_label(r'Time $t$', fontsize=plt_settings['ColorbarFont'], labelpad=10)
cbar.ax.tick_params(labelsize=plt_settings['AxesFont'])

plt.savefig(f'{plot_dir}/height_surfactant_evolution.pdf', dpi=300)
plt.close()

# Plot 2: Spacetime diagram of height
print("Creating spacetime diagram...")
times = []
height_matrix = []

# First pass: determine common x grid from first file
with open(os.path.join(output_dir, files[0])) as file:
    file.readline()
    data = np.loadtxt(file)
    x_common = np.linspace(data[:, 0].min(), data[:, 0].max(), 500)

for i, f in enumerate(files[::2]):  # Use every other file for speed
    with open(os.path.join(output_dir, f)) as file:
        header = file.readline()
        time = float(header.split('@time=')[-1])
        times.append(time)
        data = np.loadtxt(file)
        x_data = data[:, 0]
        h_data = data[:, 1]
        # Interpolate onto common grid
        h_interp = np.interp(x_common, x_data, h_data)
        height_matrix.append(h_interp)

height_matrix = np.array(height_matrix)
x = x_common

fig2, ax = plt.subplots(figsize=(12, 8))
im = ax.pcolormesh(x, times, height_matrix, shading='auto', cmap='viridis')
style_axis(ax, xlabel=r'Position $x$', ylabel=r'Time $t$', title=r'Spacetime Evolution of Film Height')
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
        gamma = data[:, 3]
        
        time_data.append(time)
        x0, h0 = find_neck(x, h)
        x0_data.append(x0)
        h0_data.append(h0)
        h_max_data.append(np.max(h))
        gamma_max_data.append(np.max(np.abs(gamma)))

fig3 = plt.figure(figsize=(14, 10))
gs = GridSpec(2, 2, figure=fig3)

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
ax4.semilogy(time_data, gamma_max_data, 'm-', linewidth=3)
style_axis(ax4, xlabel=r'Time $t$', ylabel=r'Max $|\Gamma|$',
           title=r'Maximum Surfactant Concentration')

plt.tight_layout()
plt.savefig(f'{plot_dir}/time_series_analysis.pdf', dpi=300, bbox_inches='tight')
plt.close()

# Plot 4: Zoom on the bridge region
print("Creating bridge region plot...")
fig4, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 12), sharex=True)

# Plot at different stages of coalescence
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
        gamma = data[:, 3]

        # Bridge is at x=0 where the two drops meet
        x_bridge = 0.0

        # Zoom window and sort by x-coordinate
        mask = np.abs(x - x_bridge) < 0.5
        sort_idx = np.argsort(x[mask])
        x_sorted = x[mask][sort_idx]
        h_sorted = h[mask][sort_idx]
        p_sorted = p[mask][sort_idx]
        gamma_sorted = gamma[mask][sort_idx]

        ax1.plot(x_sorted, h_sorted, color=color, linewidth=3, label=f'{label} ($t={time:.1f}$)')
        ax2.plot(x_sorted, p_sorted, color=color, linewidth=3)
        ax3.plot(x_sorted, gamma_sorted, color=color, linewidth=3)

style_axis(ax1, ylabel=r'Height $h$', title=r'Bridge Region Evolution During Coalescence')
ax1.legend(fontsize=plt_settings['LegendFont'], frameon=False)

style_axis(ax2, ylabel=r'Pressure $p$')

style_axis(ax3, xlabel=r'Position $x$', ylabel=r'Surfactant $\Gamma$')

plt.tight_layout()
plt.savefig(f'{plot_dir}/bridge_region_zoom.pdf', dpi=300, bbox_inches='tight')
plt.close()

# Plot 5: Neck trajectory (x0 vs h0)
print("Creating neck trajectory plot...")
fig5, ax = plt.subplots(figsize=(10, 10))

# Use color gradient for time
colors = plt.cm.plasma(np.linspace(0, 1, len(time_data)))

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

# Add colorbar for time
sm = plt.cm.ScalarMappable(cmap=plt.cm.plasma, norm=plt.Normalize(vmin=time_data[0], vmax=time_data[-1]))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label(r'Time $t$', fontsize=plt_settings['ColorbarFont'], labelpad=10)
cbar.ax.tick_params(labelsize=plt_settings['AxesFont'])

plt.savefig(f'{plot_dir}/neck_trajectory.pdf', dpi=300, bbox_inches='tight')
plt.close()

print(f"All plots have been saved to the '{plot_dir}' folder!")
print(f"Total number of timesteps analyzed: {len(files)}")
print(f"Time range: {time_data[0]:.2f} to {time_data[-1]:.2f}")
print(f"Final neck height: {h0_data[-1]:.6f}")
print(f"Final neck position: {x0_data[-1]:.6f}")