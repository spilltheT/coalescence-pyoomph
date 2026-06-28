# Sessile Drop Coalescence with Surfactants

Finite element simulations of viscous sessile drop coalescence in the presence of insoluble surfactants, using the lubrication approximation. This codebase accompanies a forthcoming publication in the *Journal of Fluid Mechanics*.

## Overview

This repository provides tools to simulate the early-time dynamics of two viscous sessile drops merging on a substrate. A key focus is understanding how **insoluble surfactants** modify the coalescence singularity through Marangoni stresses. The code implements a one-dimensional lubrication model coupled to interfacial surfactant transport, solved using the [pyoomph](https://pyoomph.github.io/) finite element framework.

**Key physics captured:**
- Capillary-driven bridge growth from the coalescence singularity
- Marangoni stresses arising from surfactant concentration gradients
- Competition between advective surfactant transport and surface diffusion
- Asymmetric configurations (surfactant-laden drop meets clean drop)

**Control parameters:**
- $\beta$: Surfactant strength (fractional surface tension reduction)
- $Pe$: Surface Peclet number (advection/diffusion ratio)
- $\theta$: Contact angle (sets initial drop geometry)

## Governing Equations

The dimensionless lubrication system couples film height $h(x,t)$, pressure $p(x,t)$, and surfactant concentration $\Gamma(x,t)$:

**Mass conservation:**
$$\frac{\partial h}{\partial t} + \frac{\partial}{\partial x}\left(-\frac{h^3}{3}\frac{\partial p}{\partial x} - \frac{\beta h^2}{2}\frac{\partial \Gamma}{\partial x}\right) = 0$$

**Pressure-curvature relation:**
$$p + \frac{\partial}{\partial x}\left((1 - \beta\Gamma)\frac{\partial h}{\partial x}\right) = 0$$

**Surfactant transport:**
$$\frac{\partial \Gamma}{\partial t} + \frac{\partial}{\partial x}\left(-\frac{h^2 \Gamma}{2}\frac{\partial p}{\partial x} - \beta h \Gamma \frac{\partial \Gamma}{\partial x} - \frac{1}{Pe}\frac{\partial \Gamma}{\partial x}\right) = 0$$

### Parameter Definitions

| Symbol | Name | Definition | Physical meaning |
|--------|------|------------|------------------|
| $\beta$ | Surfactant strength | $\beta = a\Gamma_\infty/\gamma_0$ | Maximum fractional reduction of surface tension |
| $Pe$ | Peclet number | $Pe = \gamma_0 L/(\mu D_s)$ | Ratio of advection to diffusion timescales |
| $\theta$ | Contact angle | Geometric parameter | Sets initial drop shape (degrees) |
| $h_\infty$ | Precursor film | Regularization | Thin film covering substrate |
| $\Gamma_0$ | Initial concentration | $\Gamma/\Gamma_\infty$ on left drop | Dimensionless surfactant loading |

### Initial Geometry

Two spherical caps touch at $x=0$:
- Drop centers at $x = \pm\sqrt{2RH - H^2}$
- Sphere radius: $R = L/\sin\theta$
- Apex height: $H = L(1-\cos\theta)/\sin\theta$
- Left drop coated with surfactant ($\Gamma = \Gamma_0$), right drop clean ($\Gamma = 0$)

## Installation

### Dependencies

- **pyoomph** (finite element framework)
- **numpy** (numerical arrays)
- **scipy** (signal processing for analysis)
- **matplotlib** (publication-quality plotting)
- **ffmpeg** (optional, for video generation)

### Installing pyoomph

**On Linux** (Python 3.9-3.13):
```bash
python -m pip install pyoomph
```
Ensure you have `gcc` installed, then verify:
```bash
python -m pyoomph check all
```

**On Mac (Apple Silicon M1-M4):**
You must run commands in a **Rosetta 2 terminal**. Also pin the MKL version:
```bash
python3 -m pip install mkl==2021.4.0
python3 -m pip install pyoomph
```
Install XCode developer tools if needed:
```bash
xcode-select --install
```

For detailed installation instructions (Windows, source compilation, troubleshooting), see the [pyoomph installation guide](https://pyoomph.github.io/installation.html).

### Setup

```bash
# Clone the repository
git clone <repository-url>
cd coalescence-with-surfactants

# Verify pyoomph installation
python -m pyoomph check all

# Test with a quick simulation
python coalescence.py --help
```

## Quick Start

Run a simulation with default parameters:

```bash
python coalescence.py
```

This creates an output directory `coalescence/` containing simulation data. Generate analysis plots:

```bash
python postprocess.py coalescence
```

Output PDFs are saved to `coalescence_plots/`.

## Simulation Scripts

### coalescence.py (Main Simulation)

Two-drop coalescence with surfactants.

**Usage:**
```bash
python coalescence.py [OPTIONS]
```

**Parameters:**

| Argument | Default | Description |
|----------|---------|-------------|
| `--beta` | 0.1 | Surfactant strength $\beta$ |
| `--Pe` | 1.0 | Peclet number $Pe$ |
| `--Gamma0` | 0.8 | Initial surfactant concentration on left drop |
| `--theta` | 20.0 | Contact angle in degrees |
| `--hp` | 1e-4 | Precursor film thickness |
| `--Lx` | 6.0 | Domain size $[-L_x/2, L_x/2]$ |
| `--N` | 1000 | Number of mesh elements |
| `--max-refinement-level` | 6 | Maximum adaptive mesh refinement level |
| `--output-dir` | coalescence | Output directory name |

**Examples:**
```bash
# Strong surfactant, low diffusion
python coalescence.py --beta 0.5 --Pe 100

# High diffusion limit
python coalescence.py --Pe 0.1 --Gamma0 0.5

# Fine mesh for convergence study
python coalescence.py --N 2000 --hp 1e-5 --output-dir coalescence_fine
```

### coalescence_clean.py

Two-drop coalescence without surfactants (baseline case).

```bash
python coalescence_clean.py
```

Uses fixed default parameters: $\theta=20°$, $h_\infty=10^{-4}$, $L_x=6$, $N=1000$.

### spreading_clean.py

Single droplet spreading on a substrate (axisymmetric).

```bash
python spreading_clean.py
```

Includes disjoining pressure for contact line dynamics.

## Postprocessing

### postprocess.py (Surfactant Case)

Generates 5 publication-quality PDF figures from surfactant simulation data.

```bash
python postprocess.py <output_directory>
python postprocess.py coalescence  # default
```

**Output plots** (saved to `<directory>_plots/`):

| File | Description |
|------|-------------|
| `height_surfactant_evolution.pdf` | Height $h(x)$ and surfactant $\Gamma(x)$ profiles at selected times |
| `spacetime_height.pdf` | Contour plot of $h(x,t)$ showing coalescence progression |
| `time_series_analysis.pdf` | 4-panel: neck height $h_0(t)$, position $x_0(t)$, $h_{\max}(t)$, $\Gamma_{\max}(t)$ |
| `bridge_region_zoom.pdf` | Detailed view of $h$, $p$, $\Gamma$ near $x=0$ at 5 stages |
| `neck_trajectory.pdf` | Phase portrait $(x_0, h_0)$ with time colorbar |

### postprocess_clean.py (Clean/Spreading Cases)

Handles both clean coalescence and axisymmetric spreading.

```bash
python postprocess_clean.py coalescence_clean
python postprocess_clean.py spreading_clean
```

**Output plots** (5-6 PDFs depending on case):
- `height_pressure_evolution.pdf`
- `spacetime_height.pdf`
- `time_series_analysis.pdf`
- `zoom_region.pdf`
- `neck_trajectory.pdf` (coalescence) or `phase_portrait.pdf` (spreading)
- `contact_line_position.pdf` (spreading only)

### postprocess_functions.py (Shared Utilities)

Core analysis functions used by all postprocessing scripts:

- `style_axis(ax, xlabel, ylabel, title)`: Publication-quality axis styling
- `find_neck(x, h)`: Robust neck detection using peak finding
- `find_drop_edge(x, h, hp, side)`: Drop boundary detection
- `plt_settings`: Font size configuration dictionary

### make_video-interfaceOnly.py

Generate animation frames from simulation data.

```bash
python make_video-interfaceOnly.py coalescence

# Encode to video with ffmpeg
ffmpeg -framerate 30 -i coalescence_frames/frame_%05d.png \
       -c:v libx264 -pix_fmt yuv420p coalescence_video.mp4
```

## Sensitivity Analysis

### check_convergence.py

Quantitative comparison of simulations with varying numerical parameters.

```bash
python check_convergence.py --hp              # Precursor film sensitivity
python check_convergence.py --Lx              # Domain size sensitivity
python check_convergence.py --N               # Mesh resolution sensitivity
python check_convergence.py --all             # All checks
python check_convergence.py --Lx --beta 0.2   # Custom parameters
```

**Output:** Formatted table with RMS% and Max% differences for $h_0(t)$, $x_0(t)$, and $x_e(t)$ relative to baseline, saved to `check-convergence-Pe{Pe}_beta{beta}-{param}.txt`.

### Shell Scripts for Parameter Sweeps

Run parallel simulations with systematic parameter variations:

```bash
./sensitivity_hp.sh              # hp = 1e-4, 4e-4, 1e-3, 1e-2
./sensitivity_Lx.sh              # Lx = 6, 8, 10, 12
./sensitivity_N.sh               # N = 500, 1000, 2000, 4000

# With custom physics parameters
./sensitivity_hp.sh --Pe 10 --beta 0.2
```

Each script:
1. Runs 4 simulations in parallel
2. Calls `check_convergence.py` for analysis
3. Generates comparison plots in `sensitivity_<param>_plots/`
4. Runs individual postprocessing for each case

## Output Format

### Directory Structure

```
coalescence/                    # Output directory
├── domain/                     # Time series data
│   ├── domain_000000.txt       # t = 0.0
│   ├── domain_000001.txt       # t = 0.1
│   └── ...
├── _ccode/                     # Generated C code (pyoomph internals)
└── _states/                    # Time-stepping state storage
```

### Text File Format

Each `domain_XXXXXX.txt` contains:

```
# coordinate_x    h    p    Gamma    @time=0.100000
-3.000000e+00    1.000000e-04    2.000000e+00    8.000000e-01
-2.994000e+00    1.000000e-04    2.000000e+00    8.000000e-01
...
```

**Columns:**
- `coordinate_x`: Spatial position
- `h`: Film height
- `p`: Pressure
- `Gamma`: Surfactant concentration (surfactant cases only)

**Header:** Contains `@time=<value>` for the simulation time.

For clean cases (no surfactant), the `Gamma` column is absent.

## Code Architecture

```
coalescence-with-surfactants/
├── Equation Modules
│   ├── lubrication.py          # LubricationEquations class (with surfactant)
│   └── lubrication_clean.py    # LubricationEquations class (clean)
│
├── Problem Definitions
│   ├── coalescence.py          # DropletCoalescence (main simulation)
│   ├── coalescence_clean.py    # DropletCoalescence (no surfactant)
│   └── spreading_clean.py      # DropletSpreading (axisymmetric)
│
├── Analysis Tools
│   ├── postprocess.py          # Surfactant case analysis
│   ├── postprocess_clean.py    # Clean/spreading analysis
│   ├── postprocess_functions.py # Shared utilities
│   ├── check_convergence.py    # Sensitivity analysis
│   └── make_video-interfaceOnly.py # Animation generation
│
└── Automation
    ├── sensitivity_hp.sh       # Precursor film sweep
    ├── sensitivity_Lx.sh       # Domain size sweep
    └── sensitivity_N.sh        # Mesh resolution sweep
```

### Class Hierarchy

- `LubricationEquations(Equations)`: Defines weak form residuals for lubrication PDEs
- `DropletCoalescence(Problem)`: Sets up mesh, initial conditions, and boundary conditions
- `DropletSpreading(Problem)`: Axisymmetric variant with disjoining pressure

## Citation

If you use this code, please cite:

```bibtex
@article{talukdar2025coalescence,
  title={Sessile-drop coalescence with surfactants},
  author={Talukdar, Jnandeep and Rocha, Duarte and Diddens, Christian
          and Snoeijer, Jacco and Lohse, Detlef and Sanjay, Vatsal},
  journal={Journal of Fluid Mechanics},
  year={2025}
}
```

## Authors

- Jnandeep Talukdar (University of Twente)
- Duarte Rocha (University of Twente)
- Christian Diddens (University of Twente)
- Jacco Snoeijer (University of Twente)
- Detlef Lohse (University of Twente / MPI Gottingen)
- Vatsal Sanjay (Durham University / CoMPhy Lab)

## License

See LICENSE file for details.
