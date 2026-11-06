# sweep

Parameter sweeps on HTCondor. Jobs save **raw** simulation output only; metrics are
derived afterwards, so they can be recomputed at any time (e.g. after a metric fix,
or for old sweeps).

```
sweep_config.yaml ──generate_params.py──► params/<run>/jobs/param_N.yaml + job_list.txt
                  ──setup_condor.sh + condor_submit──► one run_single_sim.py job per YAML
                  ──► output/sweep/<run>/jobs/<job_id>/output.pkl      (raw, one per job)
                  ──build_table.py──► output/sweep/<run>/results.csv     (one row per job)
                  ──plot_sweep.py──► output/sweep/<run>/figures/*.png
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

A copy is saved with each run (`params/<run>/` and `output/sweep/<run>/`), so a run's
data always sits next to the description of what it varied.

## generate_params.py
```
python generate_params.py <sweep_config.yaml> --name <run> [--overwrite]
```
Writes `params/<run>/jobs/param_N.yaml` (the base params with the swept keys, `SEED`,
`CB_ON` and `LOW_RES` filled in), `params/<run>/job_list.txt` (one
job YAML path per line, for Condor), and a copy of the config. Refuses an existing run
name unless `--overwrite` is passed.

## run_single_sim.py
One Condor job: simulates one job YAML and saves
`{'metadata', 'params', 'results': {t, x1, x2, spikes_n1, spikes_n2}}`. No metrics.
To test a job by hand:
```
python run_single_sim.py --params params/<run>/jobs/param_1.yaml --out-dir /tmp/test_run [--debug-plot]
```
`--debug-plot` only works locally (Condor copies back only the raw pickle).

## build_table.py
```
python build_table.py ../../output/sweep/<run> [--config PATH] [--out PATH] [--cb-on true|false]
```
Reads every `jobs/<job_id>/output.pkl`, computes metrics with
`synch_metrics.summarize_run`, and writes `results.csv`: one row per job with the swept
values (in the config's units), seed, cb_on, and chi / KOP r / firing rate /
median frequency / ISI CV for each population.

Swept values are read from what the simulation actually used: the `*_VALS` schedule
when the parameter has one, otherwise the scalar key. A schedule with more than one
entry means the parameter wasn't held fixed in that job; its value is written as NaN
with a warning.

## plot_sweep.py
```
python plot_sweep.py ../../output/sweep/<run>/results.csv [--metric r_exc r_inh ...]
                     [--sweep-params KEY ...] [--where KEY=VALUE ...] [--vmin V --vmax V]
```
Averages the realizations at each grid point and writes mean and SD figures per metric:
a line plot (1 swept parameter), a heatmap (2) or a grid of heatmaps (4). To plot two
parameters of a four-parameter sweep, fix the others, e.g.
`--sweep-params COUPLING_STRENGTH EPOP_BASE_EXCITE --where G_INTRA=0.5 --where G_INTER=1.0`.

See [../../condor/README.md](../../condor/README.md) for submitting to Condor.
