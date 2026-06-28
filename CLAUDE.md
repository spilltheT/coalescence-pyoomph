# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository implements **lubrication theory simulations for droplet coalescence with insoluble surfactants** using the [pyoomph](https://pyoomph.github.io/) finite element framework. The physics models thin liquid films on substrates where two droplets merge, optionally with surfactant transport affecting surface tension.

## Running Simulations

### Main Coalescence with Surfactants
```bash
python coalescence.py --beta 0.1 --Pe 1.0 --Gamma0 0.8 --theta 20 --hp 1e-4
```
Key parameters:
- `--beta`: Surfactant strength (fractional surface tension reduction)
- `--Pe`: Péclet number (advection/diffusion ratio)
- `--Gamma0`: Initial surfactant concentration on left droplet
- `--theta`: Contact angle in degrees
- `--hp`: Precursor film thickness
- `--Lx`: Domain size (default: 6.0, centered at x=0)
- `--N`: Number of mesh elements (default: 5000)
- `--output-dir`: Custom output directory (default: script name)

### Clean Coalescence (No Surfactant)
```bash
python coalescence_clean.py
```

### Droplet Spreading (Axisymmetric)
```bash
python spreading_clean.py
```

## Postprocessing

### Analysis Plots
```bash
python postprocess.py coalescence           # Surfactant case
python postprocess_clean.py coalescence_clean   # Clean coalescence
python postprocess_clean.py spreading_clean     # Spreading case
```

Output goes to `{base_folder}_plots/` with 5 publication-quality PDFs: height evolution, spacetime diagram, time series analysis, bridge region zoom, and neck trajectory.

### Video Generation
```bash
python make_video-interfaceOnly.py coalescence
ffmpeg -framerate 30 -i coalescence_frames/frame_%05d.png -c:v libx264 -pix_fmt yuv420p coalescence_video.mp4
```

## Architecture

### Equation Modules
- **`lubrication.py`**: Lubrication equations with surfactant transport (fields: h, p, Γ). Implements Marangoni stresses via β parameter.
- **`lubrication_clean.py`**: Clean lubrication equations without surfactant (fields: h, p). Supports optional disjoining pressure.

### Problem Classes
- **`coalescence.py`** → `DropletCoalescence`: Two spherical-cap droplets with surfactant initially on left droplet only
- **`coalescence_clean.py`** → `DropletCoalescence`: Same geometry, no surfactant
- **`spreading_clean.py`** → `DropletSpreading`: Single droplet spreading with disjoining pressure (axisymmetric)

### Shared Utilities
**`postprocess_functions.py`** provides shared analysis functions used by all postprocessing and convergence scripts:
- `style_axis(ax, xlabel, ylabel, title)`: Publication-quality axis styling (2pt spines, minor ticks, 0.3 alpha grid)
- `find_neck(x, h)`: Finds coalescence neck (local minimum closest to x=0) using `scipy.signal.find_peaks` with prominence filtering
- `find_drop_edge(x, h, hp, side)`: Finds drop edge position where h drops below threshold (2×hp or 0.01)
- `plt_settings`: Font size configuration dict for labels, axes, titles, legends, colorbars

### Physics Details
The lubrication equations (from paper §2.1-2.2):
- Mass conservation: ∂h/∂t + ∇·(h³/3·∇p + βh²/2·∇Γ) = 0
- Pressure-curvature: p = (1-βΓ)∇²h
- Surfactant transport: ∂Γ/∂t + ∇·(h²Γ/2·∇p + βhΓ∇Γ + 1/Pe·∇Γ) = 0

Initial geometry: two touching spherical caps with centers at x = ±√(2RH-H²), where R is the sphere radius derived from contact angle θ.

## Output Structure

Simulations create output directories named after the script:
```
coalescence/
├── domain/
│   ├── domain_000000.txt  # Text output: x, h, p, Gamma
│   ├── domain_000001.txt
│   └── ...
└── _ccode/                # Generated C code (pyoomph internals)
```

The text files contain columns: `coordinate_x`, `h`, `p`, `Gamma` (or just h, p for clean cases), with a header line containing `@time=<value>`.

## Sensitivity Analysis

Scripts for numerical convergence studies (4 values each, run in parallel):

```bash
./sensitivity_hp.sh   # Precursor film: hp = 1e-4, 4e-4, 1e-3, 1e-2
./sensitivity_Lx.sh   # Domain size: Lx = 6, 8, 10, 12
./sensitivity_N.sh    # Mesh resolution: N = 1000, 2000, 5000, 10000
```

Each script:
1. Runs 4 simulations in parallel
2. Calls `check_convergence.py` for quantitative analysis
3. Generates comparison plots in `sensitivity_<param>_plots/`
4. Runs individual postprocessing for each case

### Convergence Analysis Tool
```bash
python check_convergence.py --hp    # Precursor film sensitivity
python check_convergence.py --Lx    # Domain size sensitivity
python check_convergence.py --N     # Mesh resolution sensitivity
python check_convergence.py --all   # All checks
python check_convergence.py --Lx --beta 0.2 --Pe 10  # Custom parameters
```

Outputs a formatted table with RMS% and Max% differences for h₀(t), x₀(t), and xₑ(t) relative to baseline, saved to `check-convergence-Pe{Pe}_beta{beta}-{param}.txt`.

## Dependencies

- pyoomph (finite element framework)
- numpy, scipy, matplotlib (postprocessing)
- ffmpeg (video creation, optional)
