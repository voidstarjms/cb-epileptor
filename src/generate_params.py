"""
Generates params_list.txt — one line per Condor job (CE, X0 pair).
Run once before submitting: python generate_params.py
"""
import numpy as np

param1_values = np.linspace(0.0, 0.2, 8)     # CE (coupling strength)
param2_values = np.linspace(-4.5, -3.5, 8)  # X0 (epileptogenicity)

lines = []
for x0 in param2_values:
    for ce in param1_values:
        lines.append(f'{ce:.6f}, {x0:.6f}')

with open('params_list.txt', 'w') as f:
    f.write('\n'.join(lines) + '\n')

print(f"Written {len(lines)} jobs to params_list.txt")
