# ephys
Electrophysiology analysis and plotting tools.

This analyzer takes a master spreadsheet located at `cb-epileptor/ephys/master_sheet.ods`
and Igor binary files located in `cb-epileptor/ephys/binaries`. Binaries must be sorted into
subdirectories of the name format `mm dd yyyy`, and binaries must follow the name format
`R[run_num]_*_[mm]_[dd]_[yyyy].ibw` and must use underscores as their separator.

Usage: `python3 signal_processing.py -m MODE [-f FNAME] [-i IN_DIR] [--sweep SWEEP] [--sheet-path SHEET_PATH] [--show SHOW] [-v[v]] [-o OUT_DIR] [-t TYPE] [--man]`

Modes:
* single_file: Plot LFPs from a single file specified by `-f` or `--fname`.
    Pass `--sweep` to specify one sweep, otherwise all sweeps will be plotted.
* scatter_col: Plot scatter column plot of pre- and post-PEP spike counts by
    experiment type.
* mean_cv: Plot two-row bar graph of mean and CV of spike counts by experiment
    type, plotting both pre- and post-PEP.
* auto_v_man: Plot a scatter plot of the automatic spike tally against the
    manual one.
* power: Plot mean pre- and post-PEP power spectral density for a given run
    specified by `-f` or `--fname`.
* power_by_type: Plot mean pre- and post-PEP power spectral density of all runs of
    a given type specified by `-t` or `--type`.
* mean_spikes: Print mean pre- and post-PEP spike counts for each experiment type.
* wilcoxon: Print statistic value and p-value from Wilcoxon signed-rank test for
    each experiment type.
* change_score: Print pre-to-post change scores for each experiment of each type.

The default output directory for figures is `cb-epileptor/figures/ephys`
