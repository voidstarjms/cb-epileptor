# sweep

Parameter sweeps on HTCondor. Jobs save **raw** simulation output only; metrics are
derived afterwards, so they can be recomputed at any time (e.g. after a metric fix).

```
sweep_config.yaml ──generate_params.py──► data/<sweep>/jobs/param_N.yaml + job_list.txt
                  ──setup_condor.sh + condor_submit──► one run_single_sim.py job per YAML
                  ──► data/<sweep>/runs/<job_id>.pkl      (raw, one per job)
                  ──build_table.py──► data/<sweep>/results.csv    (one row per job)
                  ──plot_sweep.py──► data/<sweep>/figures/*.png
```

Everything for one sweep lives in `data/<sweep>/` at the project root:

```
data/<sweep>/
    sweep_config.yaml    copy of the config: what this run sweeps
    base_params.yaml     copy of the base params every job starts from
    job_list.txt         one job YAML path per line (relative to the project root), for Condor
    jobs/param_N.yaml    one complete params file per job
    runs/param_N.pkl     raw output, one per job
    results.csv          from build_table.py
    figures/             from plot_sweep.py
```

## sweep_config.yaml: what to sweep
The input you write. Lists the swept parameters (`sweep_params`) and run-wide settings; see
[sweep_configs/sweep_config.yaml](sweep_configs/sweep_config.yaml).

- `sweep_params`: 1 to 4 entries of `{key, unit, min, max, steps, label}`. `key` must exist in
  the base params file; `unit` is a Brian unit expression (`uS`, or `1` if unitless).
  Values are `np.linspace(min, max, steps)`. Their order is plot order: x, y, then
  panel row and panel column.
- `base_params`: params file every job starts from (relative to the config file).
- `realizations`: jobs per grid point. Realization k uses `SEED = k` at every grid point.
- `cb_on`, `low_res`: passed to `model.run_sim` in every job.

A copy is saved in each run's folder, so a run's data always sits next to the
description of what it varied.

## generate_params.py
```
python generate_params.py <sweep_config.yaml> --name <sweep> [--overwrite]
```
Creates `data/<sweep>/` with the job YAMLs (the base params with the swept keys,
`SEED`, `CB_ON` and `LOW_RES` filled in), `job_list.txt`, and copies of the config and
the base params. Refuses an existing run name unless `--overwrite` is passed, and
`--overwrite` still refuses once `runs/` holds results.

## run_single_sim.py
One Condor job: simulates one job YAML and saves
`{'metadata', 'params', 'results': {t, x1, x2, spikes_n1, spikes_n2}}`. No metrics.
To test a job by hand:
```
python run_single_sim.py --params ../../data/<sweep>/jobs/param_1.yaml --out-dir /tmp/test_run [--debug-plot]
```
`--debug-plot` only works locally (Condor copies back only the raw pickle).

## build_table.py
```
python build_table.py ../../data/<sweep> [--config PATH] [--out PATH]
```
Reads every `runs/<job_id>.pkl`, computes metrics with
`synch_metrics.summarize_run`, and writes `results.csv`: one row per job with the swept
values (in the config's units), seed, cb_on, and chi / KOP r / firing rate /
median frequency / ISI CV for each population.

## plot_sweep.py
```
python plot_sweep.py ../../data/<sweep>/results.csv [--metric r_exc r_inh ...]
                     [--vmin V --vmax V]
```
Averages the realizations at each grid point and writes mean and SD figures per metric:
a line plot (1 swept parameter), a heatmap (2) or a grid of heatmaps (4). 3 isn't supported.

See [../../condor/README.md](../../condor/README.md) for submitting to Condor.
