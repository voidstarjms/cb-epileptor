import argparse
from ephys_processing import *
import os
import sheet_parser
from scipy.stats import linregress
import sys
sys.path.append("..")
from plotting import analysis_plots

# NOTE: The default directory must exist for this to work
DEFAULT_IN_DIR = os.path.join("..", "..", "ephys", "binaries")
DEFAULT_SHEET_PATH = os.path.join("..", "..", "ephys", "master_sheet.ods")
DEFAULT_OUT_DIR = os.path.join("..", "..", "figures", "ephys")
DEFAULT_EXPERIMENT_FILE = "expt_types.csv"

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--fname', '-f', type=str, default=None,
                        help='Path to ibw binary file.')
    parser.add_argument('--in-dir', '-i', type=str, default=DEFAULT_IN_DIR,
                        help="""Path to directory of ibw binary files. 
Subdirectories must have naming scheme [m]m dd yyyy.""")
    parser.add_argument('--sweep', default=0, type=int,
                        help='The 1-indexed sweep number to examine. Leave out for all sweeps')
    parser.add_argument('--sheet-path', default=DEFAULT_SHEET_PATH, type=str,
                        help="""Path to master spreadsheet file to supply ancillary data for binary analysis.""")
    parser.add_argument('--show', default=False, action='store_true',
                        help='Show plots after completion.')
    parser.add_argument('--out-dir', '-o', default=DEFAULT_OUT_DIR, type=str,
                        help='Directory to save figures to.')
    parser.add_argument('-v', '--verbose', action='count',
                        default=0, help='Verbosity of operations. Type vv for level 2 verbosity.')
    parser.add_argument('--mode', '-m', type=str, default=None,
                        required=True,
                        help="""Plotting mode.
    \nplot_sweep: Plot LFPs from a single file specified by -f or --fname.
    Pass --sweep to specify one sweep, otherwise all sweeps will be plotted. 
    \nscatter_col: Plot scatter column plot of pre- and post-PEP spike counts by
    experiment type.
    \nmean_cv: Plot two-row bar graph of mean and CV of spike counts by experiment
    type, plotting both pre- and post-PEP.
    \nauto_v_man: Plot a scatter plot of the automatic spike tally against the
    manual one.
    \npower: Plot mean pre- and post-PEP power spectral density for a given run.
    \npower_by_type: Plot mean pre- and post-PEP power spectral density of all runs of
    a given type.
    \nmean_spikes: Print mean pre- and post-PEP spike counts for each experiment type.
    \nwilcoxon: Print statistic value and p-value from Wilcoxon signed-rank test for
    each experiment type.
    \nchange_score: Print pre-to-post change scores for each experiment of each type.
    \n""")
    parser.add_argument('--type', '-t', type=str, default=None,
                        help="Type of experiment for analysis. Used by power_by_type.")
    parser.add_argument('--man', action='store_true', default=False,
                        help="Use manual counts from master sheet.")
    parser.add_argument('--expt-list', '-e', type=str, default=DEFAULT_EXPERIMENT_FILE)

    args = parser.parse_args()
    fname = args.fname
    in_dir = args.in_dir
    sheet_path = args.sheet_path
    out_dir = args.out_dir
    sweep = args.sweep - 1 # Convert to zero-indexed
    show = args.show
    verbose = args.verbose
    mode = args.mode
    etype = args.type
    manual = args.man
    type_list_path = args.expt_list
    sheet_df = sheet_parser.parse_sheet(sheet_path)

    if mode == None:
        print("Please specify a mode with --mode.")
        sys.exit(1)

    if mode not in NO_ANALYZE_MODES:
        if mode == 'auto_v_man':
            expt_types, transients_auto = analyze_binaries(sheet_df, in_dir, type_list_path,
                                                            verbose=(verbose > 1),
                                                            aggregate=False)
            expt_types, transients_man = get_manual_transient_count(sheet_df, type_list_path,
                                                                    aggregate=False)
        else:
            if manual:
                expt_types, transients = get_manual_transient_count(sheet_df, type_list_path)
                title_suffix = "(Manual Counting)"
                out_suffix = "_man"
            else:
                expt_types, transients = analyze_binaries(sheet_df, in_dir, type_list_path,
                                                        verbose=(verbose > 1))
                title_suffix = "(Automated Counting)"
                out_suffix = "_auto"
            expt_types, split_point = eh.separate_prepost_columns(expt_types, transients)
    
    match mode:
        case 'count':
            if fname == None:
                print("Please specify a file to plot with --fname")
                sys.exit(1)
            else:
                start_sweep, idx = find_first_transient_cell_by_file(sheet_df, fname)
                prepep, postpep = get_binary_transients(fname, df=sheet_df, entry_idx=idx,
                                                        verbose=(verbose > 0),
                                                        df_start_pos=start_sweep, scan_for_nan=False)
                print("Pre-PEP transients:", prepep)
                print("Post-PEP transients:", postpep)
        case 'plot_sweep':
            if fname == None:
                print("Please specify a file to plot with --fname")
                sys.exit(1)
            else:
                get_binary_transients(fname, disp_sweep=sweep, show=show,
                                       verbose=(verbose > 0), out_dir=out_dir)
        case 'scatter_col':
            ephys_plots.plot_scatter_columns(expt_types, split_point,
                                                transients, show=show, out_dir=out_dir,
                                                out_suffix=out_suffix,
                                                title="Spikes by Experiment " + title_suffix)
        case 'mean_cv':
            ephys_plots.plot_mean_cv_bar(expt_types, split_point,
                                            transients, out_dir=out_dir,
                                            show=show, title="Spike Stats " + title_suffix,
                                            out_suffix=out_suffix)
        case 'auto_v_man':
            ephys_plots.plot_auto_v_man(transients_man, transients_auto,
                                        type_list_path=type_list_path,
                                        out_dir=out_dir, show=show)
        case 'power_raw':
            if fname == None:
                print("Please specify a file to plot with --fname")
                sys.exit(1)
            else:
                analysis_plots.plot_ephys_power_spec(out_dir, [binarywave.load(fname)['wave']['wData'][:, sweep]])
        case 'power_by_type':
            if etype == None:
                print("Please specify a type of experiment to analyze with --type or -t")
                sys.exit(1)
            else:
                prepep_array, postpep_array = get_sweeps_by_expt_type(sheet_df, in_dir, etype, verbose)
                if prepep_array.size == 0 or postpep_array.size == 0:
                    print("Specified type had 0 matches. Make sure you spelled it correctly.")
                    sys.exit(1)
                analysis_plots.plot_ephys_mean_power_spec(out_dir, [prepep_array, postpep_array], fmax=40,
                                                          fname="ephys_mean_power_"+etype,
                                                          labels=["Pre-PEP", "Post-PEP"])
        case 'power_pdf_raw':
            if fname == None:
                print("Please specify a file to plot with --fname")
                sys.exit(1)
            else:
                start_sweep_count, idx = find_first_transient_cell_by_file(sheet_df, fname)
                lfp_list = get_lfp_list(fname, df=sheet_df, entry_idx=idx,
                                        df_start_pos=start_sweep_count)
                out_name = os.path.join(out_dir, os.path.splitext(os.path.basename(fname))[0]+\
                                        "_power.pdf")
                analysis_plots.ephys_power_spec_pdf(lfp_list, fmax=100, outfile=out_name)
        case 'power_pdf':
            if fname == None:
                print("Please specify a file to plot with --fname")
                sys.exit(1)
            else:
                start_sweep_count, idx = find_first_transient_cell_by_file(sheet_df, fname)
                lfp_list = get_lfp_list(fname, df=sheet_df, entry_idx=idx,
                                        df_start_pos=start_sweep_count)
                out_name = os.path.join(out_dir, os.path.splitext(os.path.basename(fname))[0]+\
                                        "_power_normalized.pdf")
                analysis_plots.ephys_normalized_power_pdf(lfp_list, fmax=100, outfile=out_name)
        case 'power':
            if fname == None:
                print("Please specify a file to plot with --fname")
                sys.exit(1)
            else:
                start_sweep_count, idx = find_first_transient_cell_by_file(sheet_df, fname)
                lfp_list = get_lfp_list(fname, df=sheet_df, entry_idx=idx,
                                        df_start_pos=start_sweep_count)
                out_name = os.path.join(out_dir, os.path.splitext(os.path.basename(fname))[0]+\
                                        "_power_normalized.png")
                analysis_plots.ephys_normalized_power(lfp_list, fmax=30, outfile=out_name, sweep=sweep)
        case 'specgram_pdf':
            if fname == None:
                print("Please specify a file to plot with --fname")
                sys.exit(1)
            else:
                start_sweep_count, idx = find_first_transient_cell_by_file(sheet_df, fname)
                lfp_list = get_lfp_list(fname, df=sheet_df, entry_idx=idx,
                                        df_start_pos=start_sweep_count)
                out_name = os.path.join(out_dir, os.path.splitext(os.path.basename(fname))[0]+\
                                        "_sweep_spectrogram.pdf")
                ephys_plots.ephys_spectrogram_pdf(out_name, lfp_list, EPHYS_FS, fmax=10)
        case 'specgram_trace_pdf':
            if fname == None: 
                print("Please specify a file to plot with --fname")
                sys.exit(1)
            else:
                start_sweep_count, idx = find_first_transient_cell_by_file(sheet_df, fname)
                lfp_list = get_lfp_list(fname, df=sheet_df, entry_idx=idx,
                                        df_start_pos=start_sweep_count)
                out_name = os.path.join(out_dir, os.path.splitext(os.path.basename(fname))[0]+\
                                        "_sweep_spectrogram_trace.pdf")
                ephys_plots.ephys_spectrogram_trace_pdf(out_name, lfp_list, EPHYS_FS, fmax=10)
        case 'lowpass_trace':
            if sweep == -1:
                print("Please specify a sweep number with --sweep")
                sys.exit(1)
            else:
                start_sweep_count, idx = find_first_transient_cell_by_file(sheet_df, fname)
                lfp_list = get_lfp_list(fname, df=sheet_df, entry_idx=idx,
                                            df_start_pos=start_sweep_count)
                in_file_name = os.path.splitext(os.path.basename(fname))[0]
                ephys_plots.lowpass_trace(os.path.join(out_dir, f"{in_file_name}_sweep{args.sweep}_lowpass_trace.png"),
                                          lfp_list[sweep], EPHYS_FS)
        case 'esd':
            start_sweep_count, idx = find_first_transient_cell_by_file(sheet_df, fname)
            lfp_list = get_lfp_list(fname, df=sheet_df, entry_idx=idx,
                                        df_start_pos=start_sweep_count)
            ephys_plots.low_freq_esd(lfp_list, sweep, EPHYS_FS, fmax=10)
        case 'mean_spikes':
            # Print mean transient counts
            print("Mean spike counts " + title_suffix)
            print(f"{"Type":<20}{"Pre":>20}{"Post":>20}")
            for expt in expt_types:
                pre_k = expt+"_pre"
                post_k = expt+"_post"
                pre_mean = np.mean(transients[pre_k])
                post_mean = np.mean(transients[post_k])
                print(f"{expt:<20}{f"{pre_mean:.4f}":>20}{f"{post_mean:.4f}":>20}")
        case 'wilcoxon':
            print("Wilcoxon tests " + title_suffix)
            print(f"{"Type":<20}{"Statistic":>20}{"p-value":>20}")
            stat_list, pvalue_list = eh.classwise_wilcoxon(transients, expt_types)
            for expt, statistic, pvalue in zip(expt_types, stat_list, pvalue_list):
                print(f"{expt:<20}{f"{statistic:.4f}":>20}{f"{pvalue:.4f}":>20}")
        case 'change_score':
            # Print change scores for each experiment in each type
            print("Change scores " + title_suffix)
            for expt in expt_types:
                pre_k = expt+"_pre"
                post_k = expt+"_post"
                print(expt)
                for pre_ent, post_ent in zip(transients[pre_k], transients[post_k], strict=True):
                    print(f"{(post_ent+1) / (pre_ent+1) - 1:.4f} ", end='')
                print()
        case 'detect_param_sweep':
            param_step_count = 11
            detect_freq_start = 75
            detect_freq_end = 225
            detect_freq_step = (detect_freq_end - detect_freq_start) / (param_step_count - 1)
            detect_thresh_start = 6
            detect_thresh_end = 10
            detect_thresh_step = (detect_thresh_end - detect_thresh_start) / (param_step_count - 1)
            rsquare_results_pre = np.zeros((param_step_count, param_step_count))
            rsquare_results_post = np.zeros((param_step_count, param_step_count))
            freq_range = np.arange(detect_freq_start, detect_freq_end, detect_freq_step)
            thresh_range = np.arange(detect_thresh_start, detect_thresh_end, detect_thresh_step)
            _, man_transients = get_manual_transient_count(sheet_df, type_list_path, aggregate=False)
            for i, freq in enumerate(freq_range):
                for j, thresh in enumerate(thresh_range):
                    _, auto_transients = analyze_binaries(sheet_df, in_dir, type_list_path,
                                                     detect_freq=freq, detect_thresh=thresh,
                                                     verbose=(verbose > 1),
                                                     aggregate=False)
                    auto_pre, auto_post = eh.split_prepost_transients(auto_transients)
                    man_pre, man_post = eh.split_prepost_transients(man_transients)
                    _, _, r_pre, _, _ = linregress(man_pre, auto_pre)
                    _, _, r_post, _, _ = linregress(man_post, auto_post)
                    rsquare_results_pre[i, j] = r_pre**2
                    rsquare_results_post[i, j] = r_post**2
            print("Pre-PEP R squared results")
            print(rsquare_results_pre)
            max_ind = np.argmax(rsquare_results_pre)
            print("Max R squared:", np.max(rsquare_results_pre),
                  "from parameters freq", freq_range[max_ind // param_step_count],
                  "thresh", thresh_range[max_ind % param_step_count])
            print("Post-PEP R squared results")
            print(rsquare_results_post)
            max_ind = np.argmax(rsquare_results_post)
            print("Max R squared:", np.max(rsquare_results_post),
                  "from parameters freq", freq_range[max_ind // param_step_count],
                  "thresh", thresh_range[max_ind % param_step_count])
        case _:
            print("Mode not recognized. Type --help for mode list.")
               
if __name__ == "__main__":
    main()