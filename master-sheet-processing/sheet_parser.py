import pandas as pd
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--path", type=str, required=True,
                    help="Path the spreadsheet file")
args = parser.parse_args()

MASTER_SHEET_ROW_OFFSET = 23

raw_col_names = [
    "date_id_sex", "run_num", "genotype", "vehicle", "data_ok", "unused1",
    "unused2", "unused3", "total_transients", "unused4", "unused5", "unused6",
    "transient1", "transient2", "transient3", "transient4", "transient5", "transient6",
    "transient7", "transient8", "transient9", "transient10", "transient11", "transient12",
    "transient13", "transient14", "transient15", "transient16", "transient17", "transient18",
    "transient19", "transient20", "transient21", "transient22", "transient23", "transient24",
    "transient25", "transient26", "transient27", "transient28", "transient29", "transient30",
    "transient31", "transient32", "transient33", "transient34", "transient35", "transient36",
]

processed_col_names = [
    "date", "id", "sex", "run_num", "genotype", "has_cbd", "has_picro", "has_bc",
    "total_transients", "resistance", "seizure_sweep", "seizure_freq",
    "transient1", "transient2", "transient3", "transient4", "transient5", "transient6",
    "transient7", "transient8", "transient9", "transient10", "transient11", "transient12",
    "transient13", "transient14", "transient15", "transient16", "transient17", "transient18",
    "transient19", "transient20", "transient21", "transient22", "transient23", "transient24",
    "transient25", "transient26", "transient27", "transient28", "transient29", "transient30",
    "transient31", "transient32", "transient33", "transient34", "transient35", "transient36",
]

raw_dataframe = pd.read_excel(args.path, keep_default_na=False)
raw_df_view = raw_dataframe[:][MASTER_SHEET_ROW_OFFSET:]
raw_df_view.columns = raw_col_names
raw_df_view.reset_index(drop=True, inplace=True)
date_id_sex = raw_df_view["date_id_sex"][0].split()
date = date_id_sex[0]
id = date_id_sex[1]
sex = date_id_sex[2]
run_num = raw_df_view["run_num"][0]
genotype = raw_df_view["genotype"][0]
has_cbd = "CBD" in raw_df_view["vehicle"][0]
has_picro = "Picrotoxin" in raw_df_view["vehicle"][0]
has_bc = "BC" in raw_df_view["vehicle"][0]
total_transients = raw_df_view["total_transients"][2]
resistance = raw_df_view["genotype"][2]
seizure_sweep_cells = raw_df_view.iloc[0, 12:]
seizure_sweep = np.where(seizure_sweep_cells == 1)[0]
seizure_sweep = np.nan if seizure_sweep.size == 0 else seizure_sweep
if seizure_sweep is not np.nan:
    seizure_freq_cells = raw_df_view.iloc[1, 12:]
    seizure_freq = seizure_freq_cells[seizure_freq_cells.str.len() > 0]
else:
    seizure_freq = np.nan
processed_df_entry = [date, id, sex, run_num, genotype,
                            has_cbd, has_picro, has_bc, total_transients,
                            resistance, seizure_sweep, seizure_freq]
for s in range(1, 37):
    entry_name = "transient"+str(s)
    processed_df_entry.append(raw_df_view[entry_name][2])

processed_df = pd.DataFrame(columns=processed_col_names)
processed_df.loc[len(processed_df)] = processed_df_entry
print(processed_df.T)
