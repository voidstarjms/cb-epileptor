import pandas as pd
import argparse
import numpy as np
from datetime import datetime

MASTER_SHEET_ROW_OFFSET = 23
MAX_PRE_PEP_SWEEP = -17
MAX_POST_PEP_SWEEP = 24
EXPERIMENT_ENTRY_HEIGHT = 3

# Column labels for the input spreadsheet
raw_col_names = [
    "date_id_sex", "run_num", "genotype", "vehicle", "total_transients", "unused"]\
    + [("transient" + ("+" if i >= 0 else "") + str(i)) for i in
                 range(MAX_PRE_PEP_SWEEP, MAX_POST_PEP_SWEEP)]

# Column labels for the output DataFrame
processed_col_names = [
    "date", "id", "sex", "run_num", "genotype", "has_cbd", "has_picro", "has_bc",
    "data_OK", "total_transients", "resistance", "seizure_sweep", "seizure_freq"] + \
    [("transient" + ("+" if i >= 0 else "") + str(i)) for i in
     range(MAX_PRE_PEP_SWEEP, MAX_POST_PEP_SWEEP)]

def parse_sheet(path : str) -> pd.DataFrame:
    """Parse the master spreadsheet and arrange the data in a compact formatted DataFrame.
    
    Args:
        path (str): path to the spreadsheet to parse.

    Returns:
        DataFrame: the results from the master sheet formatted such that one experiment
            is one row and each column represents one field.
    """
    # Load the Excel file at the specified path into a DataFrame
    raw_dataframe = pd.read_excel(path, keep_default_na=False,
                                  na_values="#N/A")
    # Slice the raw DataFrame to only encompass the experiment data
    raw_df_view = raw_dataframe[:][MASTER_SHEET_ROW_OFFSET:]
    # Name the columns of the raw DataFrame, then reset the index
    raw_df_view.columns = raw_col_names
    raw_df_view.reset_index(drop=True, inplace=True)
    # Initialize an empty DataFrame for the processed results
    processed_df = pd.DataFrame(columns=processed_col_names)

    # Iterate rowwise down the raw DataFrame to parse the data
    ent_start = 0
    date_id_sex = []
    while ent_start < len(raw_df_view):
        # Skip blank rows (e.g. in between days or quarters)
        if raw_df_view["vehicle"][ent_start] == "":
            ent_start += 1
            continue
        
        # Update the date, ID, and sex if a new set is provided
        new_date_id_sex = raw_df_view["date_id_sex"][ent_start].split()
        if new_date_id_sex != []:
            date_id_sex = new_date_id_sex
            date_entry = datetime.strptime(date_id_sex[0], "%m/%d/%Y").date()
            id = date_id_sex[1]
            sex = date_id_sex[2]

        run_num = int(raw_df_view["run_num"][ent_start])
        genotype = raw_df_view["genotype"][ent_start]
        
        # Parse the vehicle type and set Boolean flags accordingly
        has_cbd = "cbd" in raw_df_view["vehicle"][ent_start].lower()
        has_picro = "picrotoxin" in raw_df_view["vehicle"][ent_start].lower()
        has_bc = "bc" in raw_df_view["vehicle"][ent_start].lower()

        data_ok = not raw_df_view["total_transients"][ent_start]
        total_transients = raw_df_view["total_transients"][ent_start + 2]
        resistance = raw_df_view["genotype"][ent_start + 2]

        # Record marks in the seizure marker cells
        seizure_sweep_cells = raw_df_view.iloc[ent_start, 12:]
        seizure_sweep = np.where(seizure_sweep_cells == 1)[0]
        seizure_sweep = np.nan if seizure_sweep.size == 0 else seizure_sweep
        # Record the seizure frequencies if any seizures were marked
        if seizure_sweep is not np.nan:
            seizure_freq_cells = raw_df_view.iloc[ent_start + 1, 12:]
            seizure_freq = seizure_freq_cells[seizure_freq_cells.str.len() > 0].values
        else:
            seizure_freq = np.nan
        processed_df_entry = [date_entry, id, sex, run_num, genotype,
                                    has_cbd, has_picro, has_bc, data_ok,
                                    total_transients, resistance, seizure_sweep,
                                    seizure_freq]
        
        # Fill the transient counts
        for s in range(MAX_PRE_PEP_SWEEP, MAX_POST_PEP_SWEEP):
            entry_name = "transient"+("+" if s >= 0 else "")+str(s)
            processed_df_entry.append(raw_df_view[entry_name][ent_start + 2])

        processed_df.loc[len(processed_df)] = processed_df_entry
        
        # Increment entry start position by the height of an entry (3)
        ent_start += EXPERIMENT_ENTRY_HEIGHT

    return processed_df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str, required=True,
                    help="Path the spreadsheet file")
    args = parser.parse_args()

    processed_df = parse_sheet(args.path)
    print(processed_df)

if __name__ == "__main__":
    main()