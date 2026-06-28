import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['xtick.labelsize'] = 14
matplotlib.rcParams['ytick.labelsize'] = 14

# Get base folder from command line, default to "coalescence"
base_folder = sys.argv[1] if len(sys.argv) > 1 else "coalescence"
output_dir = f"{base_folder}/domain"
frame_dir = f"{base_folder}_frames"

if not os.path.exists(frame_dir):
    os.makedirs(frame_dir)

files = sorted([f for f in os.listdir(output_dir) if f.endswith('.txt')])

print(f"Processing {len(files)} files from {output_dir}...")

for i, file in enumerate(files):
    file_path = os.path.join(output_dir, file)

    with open(file_path) as f:
        header = f.readline()
        time = float(header.split('@time=')[-1])
        data = np.loadtxt(f)

    x = data[:, 0]
    h = data[:, 1]

    # Sort by x for clean plotting
    sort_idx = np.argsort(x)
    x = x[sort_idx]
    h = h[sort_idx]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x, h, linestyle='-', color='black', linewidth=1.5)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$h$')
    ax.set_title(r'$t = %.2f$' % time)

    # Zoom around bridge region
    ax.set_xlim(-1, 1)
    ax.set_ylim(0, 0.2)
    ax.grid(True, alpha=0.3)

    output_file = os.path.join(frame_dir, f'frame_{i:05d}.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

    if (i + 1) % 100 == 0:
        print(f"  Processed {i + 1}/{len(files)} frames...")

print(f"All {len(files)} frames saved to '{frame_dir}/'")
print(f"To make a video, use:")
print(f"  ffmpeg -framerate 30 -pattern_type glob -i '{frame_dir}/*.png' -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -c:v libx264 -r 30 -pix_fmt yuv420p {base_folder}_video.mp4")
