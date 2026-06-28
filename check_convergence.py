#!/usr/bin/env python3
"""
Convergence analysis for sensitivity studies.

Usage:
    python check_convergence.py --hp    # Check precursor film sensitivity
    python check_convergence.py --Lx    # Check domain size sensitivity
    python check_convergence.py --N     # Check mesh resolution sensitivity
    python check_convergence.py --all   # Run all checks
"""

import argparse
import glob
import re
import numpy as np
import sys

from postprocess_functions import find_neck, find_drop_edge

# Display names for sensitivity parameters
PARAM_DISPLAY_NAMES = {
    'hp': 'PRECURSOR FILM THICKNESS',
    'Lx': 'DOMAIN SIZE',
    'N': 'MESH RESOLUTION',
}


def load_simulation_data(outdir, hp=1e-4):
    """
    Load h0(t), x0(t), xe(t) from simulation output directory.

    Variables:
        h0: Neck height — local minimum closest to x=0 (using find_neck)
        x0: Neck position — x-coordinate of neck minimum
        xe: Drop edge — rightmost x where h > threshold (using find_drop_edge)

    Args:
        outdir: Path to simulation output directory
        hp: Precursor film thickness (used for xe threshold)
    """
    files = sorted(glob.glob(f'{outdir}/domain/domain_*.txt'))
    if not files:
        return None

    time_data = []
    h0_data = []
    x0_data = []
    xe_data = []

    for f in files:
        # Read time from header
        with open(f, 'r') as fp:
            header = fp.readline()
            time_match = re.search(r'@time=([0-9.eE+-]+)', header)
            time = float(time_match.group(1)) if time_match else 0.0

        # Load data
        data = np.loadtxt(f, skiprows=1)
        x = data[:, 0]
        h = data[:, 1]

        # Use shared functions for robust detection
        x0, h0 = find_neck(x, h)
        xe = find_drop_edge(x, h, hp, side='right')

        time_data.append(time)
        h0_data.append(h0)
        x0_data.append(x0)
        xe_data.append(xe)

    return {
        'time': np.array(time_data),
        'h0': np.array(h0_data),
        'x0': np.array(x0_data),
        'xe': np.array(xe_data)
    }


def compute_signal_difference(t1, s1, t2, s2, n_points=200):
    """
    Compare two time signals with different time grids.

    Interpolates both signals to a common time grid and computes
    relative RMS and max differences.

    Args:
        t1, s1: Time and signal arrays for first dataset
        t2, s2: Time and signal arrays for second dataset (baseline)
        n_points: Number of points for interpolation

    Returns:
        rms_pct: Relative RMS difference (%)
        max_pct: Relative max difference (%)

    Raises:
        ValueError: If time ranges do not overlap.
    """
    # Ensure arrays are sorted by time
    idx1 = np.argsort(t1)
    t1, s1 = t1[idx1], s1[idx1]
    idx2 = np.argsort(t2)
    t2, s2 = t2[idx2], s2[idx2]

    # Find overlapping time range
    t_start = max(t1[0], t2[0])
    t_max = min(t1[-1], t2[-1])

    if t_max <= t_start:
        raise ValueError(
            f"No overlapping time range: signal 1 spans [{t1[0]:.3f}, {t1[-1]:.3f}], "
            f"signal 2 spans [{t2[0]:.3f}, {t2[-1]:.3f}]"
        )

    t_common = np.linspace(t_start, t_max, n_points)

    s1_interp = np.interp(t_common, t1, s1)
    s2_interp = np.interp(t_common, t2, s2)

    diff = s1_interp - s2_interp
    scale = np.sqrt(np.mean(s2_interp**2)) + 1e-10  # RMS of baseline

    rms_pct = np.sqrt(np.mean(diff**2)) / scale * 100
    max_pct = np.max(np.abs(diff)) / scale * 100

    return rms_pct, max_pct


def compare_datasets(data_test, data_baseline):
    """
    Compare test dataset to baseline, computing signal differences.

    Returns dict with RMS and Max differences for h0, x0, xe (all in %).
    """
    if data_test is None or data_baseline is None:
        return None

    results = {}
    for key in ['h0', 'x0', 'xe']:
        rms, maxd = compute_signal_difference(
            data_test['time'], data_test[key],
            data_baseline['time'], data_baseline[key]
        )
        results[f'{key}_rms'] = rms
        results[f'{key}_max'] = maxd

    return results


def _check_sensitivity(param_name, param_values, baseline_value, outdir_template,
                       beta, Pe, hp, hp_from_value=False):
    """
    Generic sensitivity analysis for a single parameter.

    Args:
        param_name: Name of parameter being varied (e.g., 'hp', 'Lx', 'N')
        param_values: List of parameter value strings to test
        baseline_value: Which value to use as baseline for comparison
        outdir_template: Format string for output directory (e.g., 'sensitivity_hp_{}')
        beta: Surfactant strength parameter
        Pe: Péclet number
        hp: Precursor film thickness (used for xe threshold)
        hp_from_value: If True, derive hp from each param value (for hp sensitivity)

    Returns:
        all_data: Dict mapping param values to loaded data
    """
    lines = []
    def log(msg=""):
        print(msg)
        lines.append(msg)

    display_name = PARAM_DISPLAY_NAMES.get(param_name, param_name.upper())

    log("=" * 90)
    log(f"SENSITIVITY ANALYSIS: {display_name} ({param_name})")
    if hp_from_value:
        log(f"Parameters: beta = {beta}, Pe = {Pe}")
    else:
        log(f"Parameters: beta = {beta}, Pe = {Pe}, hp = {hp}")
    log("=" * 90)

    all_data = {}

    for val in param_values:
        outdir = outdir_template.format(val)
        hp_val = float(val) if hp_from_value else hp
        log(f"\nLoading {outdir}...")
        data = load_simulation_data(outdir, hp=hp_val)
        if data is None:
            log(f"  No data found")
            continue

        all_data[val] = data
        log(f"  t_max = {data['time'][-1]:.1f}, {len(data['time'])} timesteps")

    # Compute signal differences
    if baseline_value in all_data:
        baseline = all_data[baseline_value]
        log(f"\n{'-' * 90}")
        log(f"Signal differences (baseline: {param_name} = {baseline_value})")
        log(f"{'-' * 90}")
        log(f"{param_name:<10} {'h0 RMS%':<10} {'h0 Max%':<10} {'x0 RMS%':<10} {'x0 Max%':<10} {'xe RMS%':<10} {'xe Max%':<10}")
        log(f"{'-' * 90}")

        max_h0_rms, max_h0_max = 0, 0
        max_x0_rms, max_x0_max = 0, 0
        max_xe_rms, max_xe_max = 0, 0

        for val in param_values:
            if val not in all_data or val == baseline_value:
                continue

            diff = compare_datasets(all_data[val], baseline)
            if diff is None:
                continue

            max_h0_rms = max(max_h0_rms, diff['h0_rms'])
            max_h0_max = max(max_h0_max, diff['h0_max'])
            max_x0_rms = max(max_x0_rms, diff['x0_rms'])
            max_x0_max = max(max_x0_max, diff['x0_max'])
            max_xe_rms = max(max_xe_rms, diff['xe_rms'])
            max_xe_max = max(max_xe_max, diff['xe_max'])

            log(f"{val:<10} {diff['h0_rms']:<10.2f} {diff['h0_max']:<10.2f} "
                f"{diff['x0_rms']:<10.2f} {diff['x0_max']:<10.2f} "
                f"{diff['xe_rms']:<10.2f} {diff['xe_max']:<10.2f}")

        log(f"\n{'=' * 90}")
        log(f"SUMMARY ({param_name}: {baseline_value} baseline):")
        log(f"  h0(t): max RMS = {max_h0_rms:.2f}%, max peak = {max_h0_max:.2f}%")
        log(f"  x0(t): max RMS = {max_x0_rms:.2f}%, max peak = {max_x0_max:.2f}%")
        log(f"  xe(t): max RMS = {max_xe_rms:.2f}%, max peak = {max_xe_max:.2f}%")
        log(f"{'=' * 90}")

    # Save to file
    filename = f'check-convergence-Pe{Pe}_beta{beta}-{param_name}.txt'
    with open(filename, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"\nResults saved to {filename}")

    return all_data


def check_hp_sensitivity(beta=0.1, Pe=1.0):
    """Check convergence for precursor film thickness sensitivity."""
    return _check_sensitivity(
        param_name='hp',
        param_values=['1e-4', '4e-4', '1e-3', '1e-2'],
        baseline_value='1e-4',
        outdir_template='sensitivity_hp_{}',
        beta=beta, Pe=Pe, hp=1e-4,
        hp_from_value=True
    )


def check_Lx_sensitivity(beta=0.1, Pe=1.0, hp=1e-4):
    """Check convergence for domain size sensitivity."""
    return _check_sensitivity(
        param_name='Lx',
        param_values=['6', '8', '10', '12'],
        baseline_value='6',
        outdir_template='sensitivity_Lx_{}',
        beta=beta, Pe=Pe, hp=hp
    )


def check_N_sensitivity(beta=0.1, Pe=1.0, hp=1e-4):
    """Check convergence for mesh resolution sensitivity."""
    return _check_sensitivity(
        param_name='N',
        param_values=['1000', '2000', '5000', '10000'],
        baseline_value='5000',
        outdir_template='sensitivity_N_{}',
        beta=beta, Pe=Pe, hp=hp
    )


def main():
    parser = argparse.ArgumentParser(
        description='Convergence analysis for sensitivity studies',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python check_convergence.py --hp    # Check precursor film sensitivity
    python check_convergence.py --Lx    # Check domain size sensitivity
    python check_convergence.py --N     # Check mesh resolution sensitivity
    python check_convergence.py --all   # Run all checks
    python check_convergence.py --Lx --beta 0.2 --Pe 10 --hp-val 1e-3
        """
    )
    # Sensitivity type selection
    parser.add_argument('--hp', action='store_true',
                        help='Check precursor film thickness sensitivity')
    parser.add_argument('--Lx', action='store_true',
                        help='Check domain size sensitivity')
    parser.add_argument('--N', action='store_true',
                        help='Check mesh resolution sensitivity')
    parser.add_argument('--all', action='store_true',
                        help='Run all sensitivity checks')
    # Simulation parameters
    parser.add_argument('--beta', type=float, default=0.1,
                        help='Surfactant strength (default: 0.1)')
    parser.add_argument('--Pe', type=float, default=1.0,
                        help='Péclet number (default: 1.0)')
    parser.add_argument('--hp-val', type=float, default=1e-4,
                        help='Precursor film thickness for xe threshold (default: 1e-4)')

    args = parser.parse_args()

    # Default to --all if no sensitivity type specified
    if not (args.hp or args.Lx or args.N or args.all):
        parser.print_help()
        sys.exit(1)

    if args.hp or args.all:
        check_hp_sensitivity(beta=args.beta, Pe=args.Pe)
        print("\n")

    if args.Lx or args.all:
        check_Lx_sensitivity(beta=args.beta, Pe=args.Pe, hp=args.hp_val)
        print("\n")

    if args.N or args.all:
        check_N_sensitivity(beta=args.beta, Pe=args.Pe, hp=args.hp_val)


if __name__ == '__main__':
    main()
