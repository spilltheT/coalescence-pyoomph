# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Important Guidelines

**Do not run Python scripts directly.** If script output is needed, ask the user to run the script and share the results.

## Project Overview

This repository contains analysis scripts and simulation data for studying **droplet coalescence with surfactants**. The research investigates how surfactant concentration affects the coalescence dynamics of thin liquid films, with key parameters being:
- **β (beta)**: Surfactant elasticity parameter (0.5 or 0.8)
- **Pe**: Péclet number (1, 10, etc.) - ratio of advection to diffusion
- **θ (theta)**: Contact angle (10° or 20°)
- **Γ₀ (Gamma_0)**: Initial surfactant concentration (fixed at 0.8)

## Commands

### Python Analysis
```bash
python XvsT.py              # Plot x(t) vs t (uncompensated and compensated)
python fit_prefactor.py     # Fit prefactor A for x = A·t^(3/2) scaling
```
Outputs:
- `XvsT_uncompensated.pdf`, `XvsT_compensated.pdf`
- `fit_prefactor.pdf`

**Dependencies**: numpy, pandas, matplotlib (with LaTeX), scipy

### MATLAB Analysis
```matlab
run('XvsT.m')      % Plot x(t) - coalescence position vs time
run('HvsT.m')      % Plot h(t) - film thickness vs time
run('fig9.m')      % Plot domain profiles (h, Γ vs x)
run('xfvsPe_v1.m') % Plot final position vs Péclet number
```

## Dataset Parameter Mapping

The 8 datasets in `h_f_vals/min_h_values_{1-8}.txt` correspond to specific parameter combinations:

| Dataset | β   | Pe  | θ (deg) | Time array |
|---------|-----|-----|---------|------------|
| 1       | 0.8 | 1   | 10      | t1 (0-1000, 10001 pts) |
| 2       | 0.8 | 10  | 10      | t1 |
| 3       | 0.5 | 1   | 10      | t1 |
| 4       | 0.5 | 10  | 10      | t1 |
| 5       | 0.8 | 1   | 20      | t2 (0-250, 2501 pts) |
| 6       | 0.8 | 10  | 20      | t2 |
| 7       | 0.5 | 1   | 20      | t2 |
| 8       | 0.5 | 10  | 20      | t2 |

## Data Structure

### Simulation Output (`h_f_vals/`)
- `min_h_values_{1-8}.txt`: Tab-separated with columns `File`, `Min_h`, `Min_x`
- `Min_h`: minimum film thickness at that timestep
- `Min_x`: x-position of that minimum (coalescence front position)

### Domain Snapshots (`9029/`)
- `domain_NNNNNN.txt`: Tab-separated spatial profiles at timestep NNNNNN
- Columns: `coordinate_x`, `h`, `p`, `Gamma` (header includes `@time=...`)

### Parameter Summary (`xf_vs_Pe.csv`)
- Columns: `Ma` (Marangoni number), `Pe`, `x_f` (final coalescence position)

## Physics Context

The scripts analyze power-law scaling behavior:
- **Uncompensated**: Raw x(t) or h(t) data on log-log plots
- **Compensated**: Data normalized by theoretical scaling factor `x_a` or `h_a` to collapse curves

Key theoretical scalings (Python version):
- `c_v = 0.272·(1 - β·Γ₀/2)·θ⁴`
- `x_a = c_v·β·Γ₀·√Pe / (3√π)` - normalization for position
- `h_a = (1 - β·Γ₀/2)·θ⁴` - normalization for thickness
- Expected power law: x ~ t^(3/2), h ~ t
