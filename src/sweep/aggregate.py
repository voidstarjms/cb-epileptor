"""Aggregation script — run after all Condor jobs finish.

Reads all per-job result pkls under data/results/, groups by (ce, x0, Gintra,
Ginter), computes mean and SD chi across realizations, and writes heatmap PNGs
to figures/. The run number is prefixed onto output filenames so successive
sweeps don't overwrite each other.

Usage:

    python aggregate.py             # default 2D heatmap
    python aggregate.py --full      # multifaceted (4D) panel grid
"""
import os
import sys
import argparse
import pickle
import numpy as np
from collections import defaultdict

sys.path.append("..")
from run import FilePaths
import plotting.sync_heatmap as ps

parser = argparse.ArgumentParser(description='Aggregate synchrony values from experiments into a heat map')
parser.add_argument('--full', default=False, action='store_true', help='Plot multifaceted heatmap to include synaptic conductances')
parser.add_argument('--datadir', default='data', type=str, help='Data directory to load from. Default: data')
args = parser.parse_args()
datadir = args.datadir

results_dir = os.path.join(datadir, 'results')
sidecar_path = os.path.join(datadir, 'sidecar.txt')

if not os.path.exists(results_dir):
    print(f"No results directory found at {results_dir}")
    sys.exit(1)

# Load all job result pkls
all_results = []
for fname in sorted(os.listdir(results_dir)):
    if fname.endswith('.pkl'):
        with open(os.path.join(results_dir, fname), 'rb') as f:
            contents = pickle.load(f)
            all_results.append(contents)

if len(all_results) == 0:
    print(f"No result files found in {results_dir}")
    sys.exit(1)

print(f"Loaded {len(all_results)} job results")

# Reconstruct param axes from results
J_values = sorted(set(round(r['params']['COUPLING_STRENGTH'], 6) for r in all_results))
excite_values = sorted(set(round(r['params']['EPOP_BASE_EXCITE'], 6) for r in all_results))
Gintra_values = sorted(set(round(r['params']['G_INTRA'], 6) for r in all_results))
Ginter_values = sorted(set(round(r['params']['G_INTER'], 6) for r in all_results))
print(f"J values  ({len(J_values)}): {[round(v,3) for v in J_values]}")
print(f"Epop excite values ({len(excite_values)}): {[round(v,3) for v in excite_values]}")
print(f"Gintra values ({len(Gintra_values)}): {[round(v,3) for v in Gintra_values]}")
print(f"Ginter values ({len(Ginter_values)}): {[round(v,3) for v in Ginter_values]}")

p1_label = r'Coupling Strength ($C_E$)'
p2_label = r'Innate Excitability ($x_0$)'

# Group chi values by (ce, x0, Gintra, Ginter) across realizations
chi_by_point = defaultdict(list)
for r in all_results:
    #print(r)
    key = (round(r['J'], 6), round(r['epop_excite'], 6), round(r['Gintra'], 6), round(r['Ginter'], 6))
    chi_by_point[key].append(r['r_exc'] + r['r_inh'])

figures_dir = 'figures'
os.makedirs(figures_dir, exist_ok=True)
run_num_file = 'current_run.txt'
filepaths = FilePaths(
    None, 'figures'
)

# 2-axis plot
if args.full == False:
    # Build mean and SD grids — shape: (len(x0), len(ce))
    chi_grid = np.full((len(excite_values), len(J_values)), np.nan)
    chi_sd   = np.full((len(excite_values), len(J_values)), np.nan)
    Gintra_min = min(Gintra_values)
    Ginter_min = min(Ginter_values)
    for (J, excite, Gintra, Ginter), chis in chi_by_point.items():
        # Only generate chi_grid elements for min G values for the simple graph
        if Gintra != Gintra_min or Ginter != Ginter_min:
            continue
        i = excite_values.index(excite)
        j = J_values.index(J)
        chi_grid[j, i] = np.mean(chis)
        chi_sd[j, i]   = np.std(chis)

    missing = int(np.sum(np.isnan(chi_grid)))
    if missing > 0:
        print(f"Warning: {missing} grid points are missing results")
    else:
        print("All grid points complete.")

    # Read run number set by setup_condor.sh so logs/figures/debug all share the same number
    if os.path.exists(run_num_file):
        with open(run_num_file) as f:
            run_num = int(f.read().strip())
    else:
        # Fallback: pick next available number
        run_num = 1
        while os.path.exists(os.path.join(figures_dir, f'{run_num}_sweep_debug')):
            run_num += 1


    ps.plot_synchrony(filepaths, chi_grid, chi_sd,
                  np.array(J_values), np.array(excite_values),
                  p1_label, p2_label,
                  run_num=run_num)

# 4-axis plot
else:
    p3_label = r'$G_{inter}$'
    p4_label = r'$G_{intra}$'

    # Build mean and SD grids — shape: (len(x0), len(ce), len(Gintra), len(Ginter))
    chi_grid = np.full((len(excite_values), len(J_values), len(Gintra_values), len(Ginter_values)), np.nan)
    chi_sd   = np.full((len(excite_values), len(J_values), len(Gintra_values), len(Ginter_values)), np.nan)
    for (J, excite, Gintra, Ginter), chis in chi_by_point.items():
        i = excite_values.index(excite)
        j = J_values.index(J)
        l = Gintra_values.index(Gintra)
        k = Ginter_values.index(Ginter)
        chi_grid[i, j, k, l] = np.mean(chis)
        chi_sd[i, j, k, l]   = np.std(chis)

    missing = int(np.sum(np.isnan(chi_grid)))
    if missing > 0:
        print(f"Warning: {missing} grid points are missing results")
    else:
        print("All grid points complete.")

    # Read run number set by setup_condor.sh so logs/figures/debug all share the same number
   
    if os.path.exists(run_num_file):
        with open(run_num_file) as f:
            run_num = int(f.read().strip())
    else:
        # Fallback: pick next available number
        run_num = 1
        while os.path.exists(os.path.join(figures_dir, f'{run_num}_sweep_debug_full')):
            run_num += 1

    ps.plot_synchrony_multifaceted(filepaths, chi_grid, chi_sd,
                  np.array(excite_values), np.array(J_values), np.array(Ginter_values),
                  np.array(Gintra_values), p2_label, p1_label, p3_label, p4_label)