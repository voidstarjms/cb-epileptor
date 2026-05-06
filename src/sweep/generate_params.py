"""
Run once before submitting: python generate_params.py
"""
import numpy as np

GRANULARITY = 4

param1_values = np.linspace(0.0, 0.2, GRANULARITY)    # CE (coupling strength)
param2_values = np.linspace(-4.5, -2.5, GRANULARITY)  # X0 (epileptogenicity)
param3_values = np.linspace(0.5, 2.0, GRANULARITY)      # Gintra (intrapopulation synapse weight)
param4_values = np.linspace(1.0, 4.0, GRANULARITY)        # Ginter (interpopulation synapse weight)
n_realizations = 1

lines = []
for x0 in param2_values:
    for ce in param1_values:
        for g_intra in param3_values:
            for g_inter in param4_values:
                for r in range(1, n_realizations + 1):
                    lines.append(f'{ce:.6f}, {x0:.6f}, {g_intra:.6f}, {g_inter:.6f}, {r}')

with open('params_list.txt', 'w') as f:
    f.write('\n'.join(lines) + '\n')

print(f"Written {len(lines)} jobs to params_list.txt ({len(lines)} = {len(param1_values)}x{len(param2_values)}\
        {len(param3_values)}x{len(param4_values)} grid x {n_realizations} realizations)")
