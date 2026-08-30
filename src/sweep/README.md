# sweep
HTCondor parallel computing utilities for conducting parameter sweeps.

## aggregate.py
Take in results from a given run and output a synchrony heatmap from them.

Usage: `python3 aggregate.py [--full] [--name NAME]`
if `--name` is not passed, default run name is `default`.
Passing `--full` produces a four-dimensional faceted heatmap, otherwise a
simple two-dimensional heatmap is produced.

Figures are output to `cb-epileptor/figures/sweep/[name]_synchrony_r_mean[.png | _full.png]`.

## generate_params.py
Generate a set of parameter YAML files for a Condor sweep. Parameter files make up a
Cartesian product of a set of parameters provided in the sweep manifest, swept at a
provided granularity (i.e. each parameter is interpolated from its assigned minimum
value to its maximum value over a number of steps equal to the provided granularity).

A sweep manifest is a plain text file with lines structured as
`[param_name] [axis_label]`.
Parameter names are contiguous strings (typically with underscores as separators), while
axis labels may be space-separated. In other words, the first string when the line is
split by spaces is the parameter name, while the rest of the line is the axis label.
The default sweep manifest path is `sweep_manifest.txt`.

Usage: `python3 generate_params.py -g GRANULARITY [-i INFILE] [-r REALIZATIONS] [-o OUTDIR]`.

The default output directory is `default`.

## run_single_sim.py
Run a single simulation with a given set of parameters, compute synchrony metrics on the
results, and output a small pkl file containing a summary of the results. This program is
the main driver for a single Condor job.

Usage: `python3 run_single_sim.py [--params PARAMS] [--realization REALIZATION] [--debug-plot DEBUG_PLOT]

In general, you should not have to run this program manual except for testing and debugging.

## Parameter sweep workflow (HTCondor)

One Condor job is dispatched for each YAML in `src/sweep/params/[run_name]/params`.
Each job runs a single sim, computes the chi and Kuramoto order parameter synchrony measures over
both neuron populations, and writes a compact result file. `aggregate.py` then
assembles the synchrony grid into heatmaps.

```
cb-epileptor/output/[run_name]
  data/jobs/<job_id>/output.pkl       # full per-job sim output
  data/results/<job_id>.pkl           # compact {params, chis, KOP r vals}
  figures/sweep_debug/<job_id>/       # per-job debug plots
  figures/<run>_synchrony_chi_*.png   # heatmaps from aggregate.py
```

See [condor/README.md](condor/README.md) for the full submission workflow.
