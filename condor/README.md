# HTCondor Sweep Submission

Runs a parameter sweep on the WWU cluster. See [../src/sweep/README.md](../src/sweep/README.md)
for what each sweep script does.

## Files

- `make_env_bundle.sh`: packs the venv into `cn_venv.tar.gz`, which is sent to every job.
- `setup_condor.sh`: writes `condor.sub` for one run and creates its `runs/` folder.
- `wrapper.sh`: the Condor executable. Unpacks the venv, `cd`s to `src/sweep`, runs Python.
- `logs/<sweep>/`: Condor stdout/stderr/log files.
- `condor.sub`, `cn_venv.tar.gz`: generated, not committed.

## Prerequisites

- Build `cn_venv.tar.gz` **on a pool machine**: `CN_VENV=/path/to/.cn_venv condor/make_env_bundle.sh`.
  `venv-pack` packs the installed packages but not the Python interpreter, so execute
  machines need the same base Python at the same path as the machine that built it. A
  bundle built on a Mac won't run on the Linux nodes. Rebuild whenever dependencies change.

## How to run a sweep

Commands run from the `cb-epileptor/` project root unless noted.

**1. Write a sweep config** (what to sweep): copy `src/sweep/sweep_configs/sweep_config.yaml`
and edit the sweep_params, realizations, `cb_on` and `low_res`.

**2. Generate the job YAMLs**
```
cd src/sweep
python generate_params.py path/to/sweep_config.yaml --name <sweep>
```
Creates `data/<sweep>/` (job YAMLs, `job_list.txt`, copies of the config and base params).

**3. Test one job locally**
```
python run_single_sim.py --params ../../data/<sweep>/jobs/param_1.yaml --out-dir /tmp/<sweep>_test
cd ../..
```
Confirm `/tmp/<sweep>_test/runs/param_1.pkl` appears.

**4. Write condor.sub**
```
condor/setup_condor.sh <sweep>
```
Creates `data/<sweep>/runs/` and writes `condor/condor.sub`. Refuses a run whose
`runs/` folder already holds results. Set `REQUEST_MEMORY=…` to change the per-job request
(default 32GB; check the real usage in `condor/logs/<sweep>/condor.log` and lower it).

**5. Submit**
```
condor_submit condor/condor.sub
condor_q
```
Each line of `job_list.txt` becomes one job. Condor sends each job its own job YAML,
and the job sends back its raw output to `data/<sweep>/runs/<job_id>.pkl`.

**6. Build the results table and plot**
```
cd src/sweep
python build_table.py ../../data/<sweep>
python plot_sweep.py ../../data/<sweep>/results.csv --metric r_exc r_inh
```
Figures go to `data/<sweep>/figures/`.

## Debugging a failed job

Check its error log: `condor/logs/<sweep>/err_<process>.log`. Process N runs job
`param_<N+1>` (line N+1 of `job_list.txt`). A job with no `<job_id>.pkl` in
`data/<sweep>/runs/` failed or hasn't finished; `build_table.py` skips
it and the plot shows that grid point as missing.
