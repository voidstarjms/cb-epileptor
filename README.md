# cb-epileptor

Two-population neuron model simulating epileptic behavior and endocannabinoid dynamics using Brian2.

The model couples an excitatory population (Hindmarsh-Rose model) with an
inhibitory population (Morris-Lecar model) through first-order chemical synapses.
The synapses incorporates voltage-gated bidirectional plasticity simulating
endocannabinoid retrograde signaling. The neural network reproduces seizure
dynamics, including preictal activity increase, ictal synchronization, and
subtle status epilepticus in post-ictal stage.

## Setup

```sh
python -m venv ../.cn_venv      # parent dir, matches Condor wrapper.sh expectations
source ../.cn_venv/bin/activate
pip install -r requirements.txt
```

The Condor sweep scripts hardcode the venv at `SR-CB/.cn_venv/` (one level
above this repo). If you only run sims locally you can put the venv anywhere
and activate it manually.

## Quick start

```sh
cd src
python run.py -m rp
```

Loads `params.yaml` from the repo root, runs a 120 s simulation with 10
neurons per population, and writes plots to `src/output/figures/`.

## Repository layout

```
params.yaml             # default simulation parameters (Brian2 quantity strings)
requirements.txt        # brian2, brian2tools, matplotlib, scipy, numpy, pyyaml
branching.md            # git workflow conventions

condor/                 # HTCondor sweep submission
  README.md             # full sweep workflow
  setup_condor.sh       # generates condor.sub
  wrapper.sh            # activates venv, execs python

src/
  run.py                # main entrypoint (CLI: -m rp / rpa / rpf)
  model.py              # Brian2 sim: HR + ML populations, synapses, plasticity
  param_loader.py       # YAML -> dict, resolving Brian2 expressions
  data_processing.py    # save/load output.pkl, spike histograms
  synch.py              # chi synchrony measure + Kuramoto order parameter
  plotting/             # LFP/raster, plasticity, signal-analysis, sweep heatmap plots
  sweep/                # parameter-sweep machinery (drives Condor)
    generate_params.py  # writes params/param_N.yaml + params_list.txt
    run_single_sim.py   # one Condor job: run sim, save chi summary
    aggregate.py        # collate per-job results into heatmaps
```

## Single-simulation workflow (`run.py`)

```sh
python run.py [-m MODE] [--params PATH] [--out-dir DIR] [--no-cb]
```

Mode flags are independent and stackable — e.g. `-m rpa` means run, plot,
then analyze.

| Flag | Action |
|---|---|
| `r` | Run the simulation, write `output.pkl` |
| `p` | Generate the LFP+raster plot and (if cb is on) plasticity plot |
| `f` | Combined with `p`, also draw HR (x,y,z,I_syn) and ML (x,n) traces |
| `a` | Print chi/KOP stats and write autocorr+KOP plots |

`--no-cb` keeps the Wpre dynamics evolving but drops Wpre from the
synaptic current (so plasticity has no effect on the dynamics — useful
as a control).

`--params PATH` points at any YAML matching the schema in `params.yaml`.
All keys in `run.REQUIRED_PARAMS` must be present.

## Parameter sweep workflow (HTCondor)

The sweep dispatches one Condor job per YAML in `src/sweep/params/`. Each
job runs a single sim, computes the chi synchrony measure over the HR
population, and writes a compact result file. `aggregate.py` then
assembles the chi grid into heatmaps. See
[condor/README.md](condor/README.md) for the full submission workflow.

Edit the grid (CE × X0 × Gintra × Ginter × realizations) at the top of
`src/sweep/generate_params.py` before running it.

## Configuration

`params.yaml` is a flat YAML where most values are Brian2 quantity strings
(e.g. `SIM_DURATION: "120*second"`, `G_INTER: "1*uS"`). The loader resolves
them at load time, so any valid Brian2 expression is accepted.

The four `*_VALS` arrays (`X_NAUGHT_VALS`, `COUPLING_VALS`, `G_INTER_VALS`,
`G_INTRA_VALS`) are present for time-varying schedules, but the model
equations currently read the corresponding scalar constants
(`HR_X_NAUGHT`, `COUPLING_STRENGTH`, `G_INTER`, `G_INTRA`); the sweep
overrides both for safety.

## Output

`run.py` writes:

```
out_dir/
  data/output.pkl    # {metadata, params, results} dict — see data_processing.save_data
  figures/           # PNGs (standard_plot, N1_to_1_wpre, kop, autocorr, ...)
```

The Condor sweep writes per-job artifacts under `src/sweep/`:

```
src/sweep/
  data/jobs/<job_id>/output.pkl       # full per-job sim output
  data/results/<job_id>.pkl           # compact {ce, x0, Gintra, Ginter, realization, chi}
  figures/sweep_debug/<job_id>/       # per-job debug plots
  figures/<run>_synchrony_chi_*.png   # heatmaps from aggregate.py
```

## Git

See [branching.md](branching.md) for the team's branch-and-merge conventions.

# Equations

## Neuron Populations

### Excitatory (Hindmarsh-Rose)

\begin{multline}
    {x}' = y_1 - a x^3 + b x^2 - z + I_{\text{app,1}} + J (\bar{x} - x)
    + \sigma_1 I_\text{syn,tot} + W(t)
\end{multline}
