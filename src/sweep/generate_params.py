"""
Run once before submitting: python generate_params.py
"""
import numpy as np

param1_values = np.linspace(0.0, 0.5, 8)    # CE (coupling strength)
param2_values = np.linspace(-4.5, -2.5, 8)  # X0 (epileptogenicity)
n_realizations = 5

lines = []
for x0 in param2_values:
    for ce in param1_values:
        for r in range(1, n_realizations + 1):
            lines.append(f'{ce:.6f}, {x0:.6f}, {r}')

with open('params_list.txt', 'w') as f:
    f.write('\n'.join(lines) + '\n')

print(f"Written {len(lines)} jobs to params_list.txt ({len(lines)} = {len(param1_values)}x{len(param2_values)} grid x {n_realizations} realizations)")
