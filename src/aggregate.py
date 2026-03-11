"""
Aggregation script — run after all Condor jobs finish.
Reads all per-job result pkls, reconstructs chi grid, generates heatmap plots.
Run as: python aggregate.py
"""
import os
import sys
import pickle
import numpy as np
from collections import defaultdict

import plot_synchrony as ps

results_dir = os.path.join('data', 'results')

if not os.path.exists(results_dir):
    print(f"No results directory found at {results_dir}")
    sys.exit(1)

# Load all job result pkls
all_results = []
for fname in sorted(os.listdir(results_dir)):
    if fname.endswith('.pkl'):
        with open(os.path.join(results_dir, fname), 'rb') as f:
            all_results.append(pickle.load(f))

if len(all_results) == 0:
    print("No result files found in data/results/")
    sys.exit(1)

print(f"Loaded {len(all_results)} job results")

# Reconstruct param axes from results
ce_values = sorted(set(round(r['ce'], 6) for r in all_results))
x0_values = sorted(set(round(r['x0'], 6) for r in all_results))
print(f"CE values  ({len(ce_values)}): {[round(v,3) for v in ce_values]}")
print(f"X0 values ({len(x0_values)}): {[round(v,3) for v in x0_values]}")

# Group chi values by (ce, x0) across realizations
chi_by_point = defaultdict(list)
for r in all_results:
    key = (round(r['ce'], 6), round(r['x0'], 6))
    chi_by_point[key].append(r['chi'])

# Build mean and SD grids — shape: (len(x0), len(ce))
chi_grid = np.full((len(x0_values), len(ce_values)), np.nan)
chi_sd   = np.full((len(x0_values), len(ce_values)), np.nan)
for (ce, x0), chis in chi_by_point.items():
    i = ce_values.index(ce)
    j = x0_values.index(x0)
    chi_grid[j, i] = np.mean(chis)
    chi_sd[j, i]   = np.std(chis)

missing = int(np.sum(np.isnan(chi_grid)))
if missing > 0:
    print(f"Warning: {missing} grid points are missing results")
else:
    print("All grid points complete.")

# Read run number set by setup_condor.sh so logs/figures/debug all share the same number
figures_dir = 'figures'
os.makedirs(figures_dir, exist_ok=True)
run_num_file = 'current_run.txt'
if os.path.exists(run_num_file):
    with open(run_num_file) as f:
        run_num = int(f.read().strip())
else:
    # Fallback: pick next available number
    run_num = 1
    while os.path.exists(os.path.join(figures_dir, f'{run_num}_sweep_debug')):
        run_num += 1

p1_label = r'$C_E$ (coupling strength)'
p2_label = r'$x_0$ (epileptogenicity)'

ps.plot_synchrony(chi_grid, chi_sd,
                  np.array(ce_values), np.array(x0_values),
                  p1_label, p2_label,
                  run_num=run_num)

print(f"Saved: {run_num}_synchrony_chi_mean.png and {run_num}_synchrony_chi_sd.png")
