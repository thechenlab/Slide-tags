import os
import re
import csv
import sys
import warnings
import pandas as pd
from googleapiclient.discovery import build
from google.oauth2.service_account import Credentials

warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None


sheet_id = sys.argv[1]
bcl = sys.argv[2]
src_path = sys.argv[3]
data_path = sys.argv[4]
default_bcl_path = sys.argv[5]
reference_path = sys.argv[6]


# ------------------ Parse the google sheet ------------------
def pars_google_sheet(src_path):
    """Parse the google sheet"""
    # Get the google sheet
    workflow_path = src_path.replace("scripts", "config")
    service_account_json_file = os.path.join(workflow_path, 'google_key.json')
    SERVICE_ACCOUNT_FILE = service_account_json_file
    SCOPES = ['https://www.googleapis.com/auth/spreadsheets']
    creds = Credentials.from_service_account_file(SERVICE_ACCOUNT_FILE, scopes=SCOPES)
    service = build('sheets', 'v4', credentials=creds)
    return service


# ------------------ Match BCL in spreadsheet ------------------
# Match BCL in spreadsheet
def match_bcl_arg(bcl, service, sheet_id):
    """Match BCL in spreadsheet"""
    found = False 
    SAMPLE_SPREADSHEET_ID = sheet_id
    for tech in ["Tags", "Part"]:
        SAMPLE_RANGE_NAME = tech
        sheet = service.spreadsheets()
        result = sheet.values().get(spreadsheetId=SAMPLE_SPREADSHEET_ID, range=SAMPLE_RANGE_NAME).execute()
        values = result.get('values', [])
        if not values:
            continue  
        df = pd.DataFrame(values, columns=values[0]).iloc[1:]
        rows = df[(df['Run'] == 'YES') & (df['BCL'] == bcl)]
        if len(rows) > 0:
            found = True  
            break  
    if not found:
        raise ValueError(f"ERROR: input BCL '{bcl}' had no spreadsheet matches in Tags. Please check the BCL name and try again.")
    return rows, df[df['BCL'] == bcl]


# ------------------ Add Lane info ------------------
# Add multilanes if one sample came from different lanes
def add_multilanes(row, type, reshaped_data):
    if type == 'RNA':
        sample = row['Index_ND']
        index = row['RNA Index']
        lane = row['Lane']
    elif type == 'SB':
        sample = row['SB_ND']
        index = row['SB Index']
        lane = row['SB Lane']
    elif type == 'ATAC':
        sample = row['ATAC Name']
        index = row['ATAC Index']
        lane = row['ATAC Lane']
    if "," in lane:
        lane_subs = [int(i) for i in lane.split(",")]
        for lane in lane_subs:
            entry = {'Lane': lane, 'Sample': sample, 'Index': index}
            if entry not in reshaped_data:
                reshaped_data.append(entry)
    else:
        lane = int(lane)
        entry = {'Lane': lane, 'Sample': sample, 'Index': index}
        if entry not in reshaped_data:
            reshaped_data.append(entry)


# ------------------ make input input for bcl2fastq ------------------
def make_bcl2fastq_input(input_index, lane_number, type, bcl2fastq_folder, reference_path):
    index_info_folder = os.path.join(reference_path, 'Index_Info')
    index_files = os.listdir(index_info_folder)
    atac_index_files = [f for f in index_files if f.startswith('Single')]
    rna_index_files = [f for f in index_files if f.startswith('Dual')]
    sub_folder = os.path.join(bcl2fastq_folder, 'split_lanes')
    os.makedirs(sub_folder, exist_ok=True)
    if type == 'ATAC':
        index_files_search = atac_index_files
        file_types = 'Single'
        split_name = 'ATAC_lane'
        header_names = ["Lane,Sample_ID,index"]
        sheet_name = os.path.join(bcl2fastq_folder, "ATAC_single_index.csv")
    else:
        index_files_search = rna_index_files
        file_types = 'Dual'
        split_name = 'Index_lane'
        header_names = ["Lane,Sample_ID,index,index2"]
        sheet_name = os.path.join(bcl2fastq_folder, "RNA_dual_index.csv")
    matched_index = []
    missing_index = []
    for i, lane in zip(input_index, lane_number):
        index_found = False
        for file in index_files_search:
            index_df = pd.read_csv(os.path.join(index_info_folder, file))
            matching_rows = index_df[index_df['index_name'] == i]
            if not matching_rows.empty:
                for row in matching_rows.values:
                    if file_types == 'Dual':
                        matched_index.append([lane] + [list(row)[0], list(row)[1], list(row)[3]])
                    elif file_types == 'Single':
                        matched_index.append([lane] + [list(row)[0], list(row)[1]]) 
                        matched_index.append([lane] + [list(row)[0], list(row)[2]]) 
                        matched_index.append([lane] + [list(row)[0], list(row)[3]]) 
                        matched_index.append([lane] + [list(row)[0], list(row)[4]]) 
                index_found = True
        if not index_found:
            missing_index.append(i)
    if matched_index:
        # write all index info
        all_matched_rows = pd.DataFrame(matched_index, columns=header_names[0].split(',')).drop_duplicates().reset_index(drop=True)
        if all_matched_rows.shape[0] == 0:
            raise ValueError(f"ERROR: No index found for {type} in the index files")
        matched_index_str = all_matched_rows.to_csv(index=False, header=False).strip().split('\n')
        content = ["[Data]"] + header_names + matched_index_str
        with open(sheet_name, 'w') as file:
            file.write('\n'.join(content) + '\n')
        # split by lane
        for lane, group in all_matched_rows.groupby('Lane'):
            sub_sheet_name = os.path.join(sub_folder, f"{split_name}{lane}.csv")
            content = ["[Data]"] + header_names + group.to_csv(index=False, header=False).strip().split('\n')
            with open(sub_sheet_name, 'w') as file:
                file.write('\n'.join(content) + '\n')


# ------------------ Check duplicates ------------------
def find_duplicates(input_list):
    seen = set()
    duplicates = set()
    for item in input_list:
        if item in seen:
            duplicates.add(item)
        else:
            seen.add(item)
    return list(duplicates)


# ------------------ Save RNAcounts sublist ------------------
def save_sub_tables(row, base_name, path, col):
    if len(row) > 0:
        # get cellranger version
        if row['Cellranger']:
            cellranger = row['Cellranger']
        else:
            cellranger = 'V7'
        if row['Chemistry'] in ['5P-PE-v3', '5P-R2-v3']:
            cellranger = 'V8'
        # save the file
        if len(row[col]) != 0:
            save_name = f"{base_name}_{row[col]}.csv"
            full_path = os.path.join(path, save_name)
            with open(full_path, 'w') as file:
                file.write(f"{row[col]},{row['Transcriptome']},{row['Chemistry']},{row['VDJ']},{cellranger}\n")



# ------------------ Check file exists and get paras ------------------
def get_info_and_create_files(rows, df, bcl, data_path, default_bcl_path, reference_path):
    # pre defined
    empty = ['X', 'x', ' ', '', None, float('nan')]
    default_path = default_bcl_path
    default_path_main = data_path
    base_log_path = os.path.join(default_path_main, bcl, 'log')

    with open(os.path.join(base_log_path, 'working') , 'r') as file:
        folder_to_create = os.path.join(file.readline().strip(), 'input') 
    sys.stdout.write(f"Input sheets saved in: {folder_to_create}\n")
    
    # Check main BCL path and files
    bcl_path = rows['BCL Path'].str.strip().to_list()
    if len(set(bcl_path)) > 1:
        raise ValueError(f"ERROR: multiple BCL paths found in same main BCL: {bcl_path}")
    else: bcl_path = bcl_path[0]

    if bcl_path.startswith('/'):
        if os.path.exists(bcl_path):
            bcl_path = bcl_path
        elif os.path.exists(os.path.join(default_path, bcl_path)):
            bcl_path = os.path.join(default_path, bcl_path)
    else:
        if os.path.exists(os.path.join(default_path, bcl_path)):
            bcl_path = os.path.join(default_path, bcl_path)
    if bcl_path in empty:
        sys.stdout.write(f"WARNING: BCL path: '{bcl_path}' does not exist.")
    else:
        sys.stdout.write(f"BCL Path: {bcl_path}\n")

    # Store Path
    store_paths = rows['Store Path'].str.strip().tolist()
    for i, path in enumerate(store_paths):
        if path in empty:
            store_paths[i] = None
        elif not os.path.exists(path):
            new_path = os.path.join(default_path, path)
            if os.path.exists(new_path):
                store_paths[i] = new_path
            else:
                raise ValueError(f"ERROR: Store Path: '{new_path}' does not exist")
    rows.loc[:, 'Store Path'] = store_paths

    # Custom Sample Name
    name = rows['Sample Name'].str.strip().tolist()
    if any(n in [' ', '', None] for n in name):
        raise ValueError("ERROR: Some Sample Names are empty")
    duplicates = find_duplicates(name)
    if duplicates:
        raise ValueError(f"ERROR: Duplicate Sample Name found in spreadsheet: {', '.join(duplicates)}")

    # Check RNA Index, Lane, Species and VDJ
    for index, row in rows.iterrows():
        rna_index_empty = row['RNA Index'] in empty
        lane_empty = row['Lane'] in empty
        atac_index_empty = row['ATAC Index'] in empty
        atac_lane_empty = row['ATAC Lane'] in empty
        species_empty = row['Species'] in empty
        vdj_chain = row['VDJ']
        sb_index_empty = row['SB Index'] in empty
        sb_lane_empty = row['SB Lane'] in empty
        puck_empty = row['Puck ID'] in empty
        add_rna_index_empty = row['Add RNA Index'] in empty
        add_sb_index_empty = row['Add SB Index'] in empty
        add_puck_id_empty = row['Add Puck ID'] in empty
        add_atac_index_empty = row['Add ATAC Index'] in empty
        
        # RNA Index, Lane, Species, VDJ
        if (rna_index_empty and not lane_empty) or (not rna_index_empty and lane_empty):
            raise ValueError(f"ERROR: Inconsistency in 'RNA Index' and 'Lane' in sample '{row['Sample Name']}'")
        if (not rna_index_empty) and (species_empty or lane_empty):
            raise ValueError(f"ERROR: Missing 'Species' or 'Lane' while 'RNA Index' is present in sample '{row['Sample Name']}'")
        if vdj_chain:
            assert (not rna_index_empty) and (not lane_empty) and (not species_empty), f"ERROR: 'VDJ' is present but 'RNA Index', 'Lane', or 'Species' is missing in sample '{row['Sample Name']}'"
        
        # additional RNA index
        if (not add_rna_index_empty) and species_empty:
            raise ValueError(f"ERROR: Missing 'Species' while 'Add RNA Index' is present in sample '{row['Sample Name']}'")
        
        # additional SB index
        if (not add_sb_index_empty) and add_puck_id_empty:
            raise ValueError(f"ERROR: Missing 'Add Puck ID' while 'Add SB Index' is present in sample '{row['Sample Name']}'")
        if add_sb_index_empty and (not add_puck_id_empty):
            raise ValueError(f"ERROR: Missing 'Add SB Index' while 'Add Puck ID' is present in sample '{row['Sample Name']}'")

        # check SB Index, Lane, Puck ID
        should_raise_error = False
        if (not sb_index_empty) and (not sb_lane_empty) and (not puck_empty):
            pass 
        elif sb_index_empty and sb_lane_empty and puck_empty:
            pass  
        else:
            should_raise_error = True
            error_message = f"ERROR: 'SB Index', 'SB Lane', and 'Puck ID' must either all be empty or all be filled for sample '{row['Sample Name']}'"
        if should_raise_error:
            raise ValueError(error_message)
        
        # check ATAC Index, Lane
        should_raise_error = False
        if (not atac_lane_empty) and (not atac_index_empty):
            pass  
        elif atac_lane_empty and atac_index_empty:
            pass 
        else:
            should_raise_error = True
            error_message = f"ERROR: 'ATAC Index' and 'ATAC Lane' must either all be empty or all be filled for sample '{row['Sample Name']}'"
        if should_raise_error:
            raise ValueError(error_message)
        
        # additional ATAC index
        if (not add_atac_index_empty) and (not atac_index_empty):
            raise ValueError(f"ERROR: 'Add ATAC Index' and 'ATAC Index' should not fill at same time '{row['Sample Name']}'")


    # Add additional RNA and SB Index if not present in main BCL
    df.loc[df['RNA Index'].isin(empty) & ~df['Add RNA Index'].isin(empty), 'RNA Index'] = df['Add RNA Index']
    df.loc[df['SB Index'].isin(empty) & ~df['Add SB Index'].isin(empty), 'SB Index'] = df['Add SB Index']

    # Make unique Index within same BCL
    # RNA
    df['Index_ND'] = df['RNA Index'].copy()
    not_in_empty = ~df['Index_ND'].isin(empty)
    duplicates = df[not_in_empty].duplicated(subset='Index_ND', keep=False)
    duplicates_index = df[not_in_empty][duplicates].index
    df.loc[duplicates_index, 'Index_ND'] += df.loc[duplicates_index, 'Lane'].astype(str).apply(lambda x: '' if x == '' or ',' in x else '-' + x)
    # SB
    df['SB_ND'] = df['SB Index'].copy()
    not_in_empty = ~df['SB_ND'].isin(empty)
    duplicates = df[not_in_empty].duplicated(subset='SB_ND', keep=False)
    duplicates_index = df[not_in_empty][duplicates].index
    df.loc[duplicates_index, 'SB_ND'] += df.loc[duplicates_index, 'SB Lane'].astype(str).apply(lambda x: '' if x == '' or ',' in x else '-' + x)
    # ATAC
    df['ATAC Name'] = df['ATAC Index'].copy()
    not_in_empty = ~df['ATAC Name'].isin(empty)
    duplicates = df[not_in_empty].duplicated(subset='ATAC Name', keep=False)
    duplicates_index = df[not_in_empty][duplicates].index
    df.loc[duplicates_index, 'ATAC Name'] += df.loc[duplicates_index, 'ATAC Lane'].astype(str).apply(lambda x: '' if x == '' or ',' in x else '-' + x)

    # susbet
    df = df[df['Run'] == 'YES']
    rows.loc[:, 'Index_ND'] = df.loc[:, 'Index_ND']
    rows.loc[:, 'SB_ND'] = df.loc[:, 'SB_ND']
    rows.loc[:, 'ATAC Name'] = df.loc[:, 'ATAC Name']
    rows.loc[~df['Add SB Index'].isin(empty), 'Puck ID'] = rows.loc[:, 'Add Puck ID']
    sys.stdout.write(f"RNA Indexes: {', '.join(set([index for index in rows['Index_ND'].str.strip().dropna() if index]))}\n")
    sys.stdout.write(f"SB Indexes: {', '.join(set([index for index in rows['SB_ND'].str.strip().dropna() if index]))}\n")
    sys.stdout.write(f"Puck ID: {', '.join(set([index for index in rows['Puck ID'].str.strip().dropna() if index]))}\n")
    
    # Add Transcriptome
    rows['Transcriptome'] = None
    for index, row in rows.iterrows():
        species_empty = row['Species']
        vdj_chain = row['VDJ']
        if species_empty == 'Human':
            if vdj_chain:
                rows.loc[index, 'Transcriptome'] = 'refdata-cellranger-vdj-GRCh38'
            elif row['ATAC Index'] or row['Add ATAC Index'] not in empty:
                rows.loc[index, 'Transcriptome'] = 'refdata-arc-GRCh38-2020-A'
            else:
                rows.loc[index, 'Transcriptome'] = 'refdata-gex-GRCh38-2024-A'
        elif species_empty == 'Mouse':
            if vdj_chain:
                rows.loc[index, 'Transcriptome'] = 'refdata-cellranger-vdj-GRCm38'
            elif row['ATAC Index'] or row['Add ATAC Index'] not in empty:
                rows.loc[index, 'Transcriptome'] = 'refdata-arc-mm10-2020-A'
            else:
                rows.loc[index, 'Transcriptome'] = 'refdata-gex-mm10-2020-A'

    # Get Contact Email
    email = list(set([item for item in rows['Email'].tolist() if item not in empty]))
    if email:
        sys.stdout.write(f"Email: {', '.join(email)}\n")

    # create folder
    def check_directory(directory):
        if not os.path.exists(directory):
            os.makedirs(directory)

    sub_count_folder = os.path.join(folder_to_create, 'split_counts')
    sub_mkfastq_folder = os.path.join(folder_to_create, 'split_mkfastq')
    sub_spatial_folder = os.path.join(folder_to_create, 'split_spatial')
    bcl2fastq_folder = os.path.join(folder_to_create, 'bcl2fastq')
    check_directory(folder_to_create)
    check_directory(sub_count_folder)
    check_directory(sub_mkfastq_folder)
    check_directory(sub_spatial_folder)
    check_directory(bcl2fastq_folder)

    # RNA index
    indexes_csv_path = os.path.join(folder_to_create, "Indexes.csv")
    # Multiome ATAC index
    atac_csv_path = os.path.join(folder_to_create, "ATAC_Index.csv")
    # Single ATAC index
    atac_single_csv_path = os.path.join(folder_to_create, "ATAC.csv")
    # cellranger and cellbender inputs
    rna_counts_csv_path = os.path.join(folder_to_create, "RNAcounts.tsv")
    # SBcounts and Spatial
    sb_counts_csv_path = os.path.join(folder_to_create, "SBcounts.tsv")
    spatial_csv_path = os.path.join(folder_to_create, "Spatial.tsv")
    name_to_index_csv_path = os.path.join(folder_to_create, "name_to_index.csv")
    # additional merge files
    merge_csv_path = os.path.join(folder_to_create, "merge.csv")
    custom_refdata_csv_path = os.path.join(folder_to_create, "custom_refdata.csv")
    
    # Write the Indexes.csv file for RNA and SB
    # for cellranger mkfastq
    reshaped_data = []
    for index, row in rows.iterrows():
        if (row['RNA Index'] not in empty) and (row['Lane'] not in empty):
            add_multilanes(row, 'RNA', reshaped_data)
        if (row['SB Index'] not in empty) and (row['SB Lane'] not in empty):
            add_multilanes(row, 'SB', reshaped_data)
    dt = pd.DataFrame(reshaped_data)
    if len(dt.index) == 0:
        sys.stdout.write("Not enough information to create RNA Indexes.csv, writing a blank file... \n")
        open(indexes_csv_path, 'w').close()
    else:
        for lane, group in dt.groupby('Lane'):
            group.to_csv(os.path.join(sub_mkfastq_folder, f'Index_lane{lane}.csv'), index=False)
        dt.to_csv(indexes_csv_path, index=False)
        sys.stdout.write(f"Indexes.csv ({len(dt.index)} lines written) \n")
        # for bcl2fastq
        input_index = dt['Index'].tolist()
        lane_number = dt['Lane'].tolist()
        make_bcl2fastq_input(input_index, lane_number, 'RNA', bcl2fastq_folder, reference_path)

    # Write the ATAC_Indexes.csv file for ATAC
    reshaped_data = []
    reshaped_data_single = []
    for index, row in rows.iterrows():
        if (row['ATAC Index'] not in empty) and (row['ATAC Lane'] not in empty):
            if (row['RNA Index'] in empty) and (row['Lane'] in empty):
                add_multilanes(row, 'ATAC', reshaped_data_single)
            else:
                add_multilanes(row, 'ATAC', reshaped_data)
    dt = pd.DataFrame(reshaped_data)
    if len(dt.index) > 0:
        dt.to_csv(atac_csv_path, index=False)
        for lane, group in dt.groupby('Lane'):
            group.to_csv(os.path.join(sub_mkfastq_folder, f'ATAC_Index_lane{lane}.csv'), index=False)
        sys.stdout.write(f"ATAC_Index.csv ({len(dt.index)} lines written) \n")
        input_index = dt['Index'].tolist()
        lane_number = dt['Lane'].tolist()
        make_bcl2fastq_input(input_index, lane_number, 'ATAC', bcl2fastq_folder, reference_path)
    dt = pd.DataFrame(reshaped_data_single)
    if len(dt.index) > 0:
        dt.to_csv(atac_single_csv_path, index=False)
        for lane, group in dt.groupby('Lane'):
            group.to_csv(os.path.join(sub_mkfastq_folder, f'ATAC_lane{lane}.csv'), index=False)
        sys.stdout.write(f"ATAC.csv (ATAC only) ({len(dt.index)} lines written) \n")
        input_index = dt['Index'].tolist()
        lane_number = dt['Lane'].tolist()
        make_bcl2fastq_input(input_index, lane_number, 'ATAC', bcl2fastq_folder, reference_path)

    # Write libraries for multiome
    df = rows[rows['Run'] == 'YES'].copy()
    for index, row in df.iterrows():
        bcl1 = row['BCL']
        bcl2 = row['BCL']
        if not(row['Add ATAC Index'] in empty):
            row['ATAC Name'] = row['Add ATAC Index']
            bcl2 = row['Merge ATAC From BCL']
        if not(row['Add RNA Index'] in empty):
            row['Index_ND'] = row['Add RNA Index']
            bcl1 = row['Merge RNA From BCL']
        atac_name = row['ATAC Name']
        index_nd = row['Index_ND']
        if not(atac_name in empty) and not(index_nd in empty):
            file_name = os.path.join(sub_count_folder, f"multiome_{index_nd}.csv")
            content = [
                "fastqs,sample,library_type",
                f"{os.path.join(default_path_main, bcl1, 'fastq', str(index_nd))},{index_nd},Gene Expression",
                f"{os.path.join(default_path_main, bcl2, 'fastq', atac_name)},{atac_name},Chromatin Accessibility"
            ]
            with open(file_name, 'w') as f:
                for line in content:
                    f.write(line + '\n')
            sys.stdout.write(f"multiome libraries ({len(df.index)} lines written) \n")
    
    # Write the RNAcounts.tsv file
    if len(set([i for i in rows['Index_ND'].tolist() if i not in empty])) == 0 and len(set([i for i in rows['ATAC Index'].tolist() if i not in empty])) == 0:
        sys.stdout.write(f"Not enough information to create RNAcounts.tsv, writing a blank file... \n")
        open(rna_counts_csv_path, 'w').close()
    else:
        reshaped_data = []
        for index, row in rows.iterrows():
            if row['Index_ND'] not in empty and row['Transcriptome'] not in empty:
                reshaped_data.append({'Sample': row['Index_ND'], 'Transcriptome': row['Transcriptome']})
                if row["VDJ"] in empty:
                    save_sub_tables(row, "count", sub_count_folder, 'Index_ND')  
                else:
                    save_sub_tables(row, "vdj", sub_count_folder, 'Index_ND')
            if row['ATAC Index'] not in empty and row['Index_ND'] in empty:
                reshaped_data.append({'Sample': row['ATAC Index'], 'Transcriptome': row['Transcriptome']})
                save_sub_tables(row, "count", sub_count_folder, 'ATAC Index')
        df = pd.DataFrame(reshaped_data)
        df.to_csv(rna_counts_csv_path, sep='\t', index=False, header=False)
        sys.stdout.write(f"RNAcounts.tsv ({len(df.index)} lines written) \n")
    
    # Write the SBcounts.tsv file
    if len(set([i for i in rows['SB_ND'].tolist() if i not in empty])) == 0:
        sys.stdout.write("Not enough information to create SBcounts.tsv and Spatial.tsv, writing a blank file... \n")
        open(sb_counts_csv_path, 'w').close()
        open(spatial_csv_path, 'w').close()
    else:
        reshaped_data_sb_counts = []
        reshaped_data_spatial = []
        for index, row in rows.iterrows():
            if row['SB_ND'] not in empty:
                reshaped_data_sb_counts.append({'Sample': row['SB_ND']})
                reshaped_data_spatial.append({'Sample': row['SB_ND'], 'Puck': row['Puck ID']})
        # SBcounts.tsv
        df_sb_counts = pd.DataFrame(reshaped_data_sb_counts)
        df_sb_counts.to_csv(sb_counts_csv_path, sep='\t', index=False, header=False)
        sys.stdout.write(f"SBcounts.tsv ({len(df_sb_counts.index)} lines written) \n")

        # Split spatial SBcounts tsv
        for index, row in df_sb_counts.iterrows():
            row.to_csv(os.path.join(sub_spatial_folder, f'SBcount_{row[0]}.csv'), index=False, header=False)

        # Spatial.tsv
        df_spatial = pd.DataFrame(reshaped_data_spatial)
        df_spatial.to_csv(spatial_csv_path, sep='\t', index=False, header=False)
        sys.stdout.write(f"Spatial.tsv ({len(df_spatial.index)} lines written) \n")

    # write name to index
    df = rows[['Sample Name', 'Index_ND', 'VDJ', 'SB_ND', 'Puck ID', 'Store Path']]
    df = df.rename(columns={'Sample Name': 'Name', 'Index_ND':'RNA Index', 'VDJ':'VDJ Chain', 'SB_ND':'SB Index'})
    df = df[~df['RNA Index'].isin(empty)]
    df.replace(empty, 'X', inplace=True)
    df.to_csv(name_to_index_csv_path, index=False)

    # write merge.csv file 
    df = rows[rows.apply(lambda row: row['Merge RNA From BCL'] not in empty or row['Merge Spatial From BCL'] not in empty, axis=1)]
    df = df[['BCL', 'Merge RNA From BCL', 'Add RNA Index', 'Merge Spatial From BCL', 'Add SB Index', 'Add Puck ID', 'Sample Name', 'Index_ND', 'SB_ND']]
    if len(df.index) > 0:
        df.to_csv(merge_csv_path, index=False)
        sys.stdout.write(f"merge.csv ({len(df.index)} lines written) \n")

    # write custom reference file
    df = rows[rows.apply(lambda row: row['Custom Reference'] not in empty, axis=1)]
    if len(df.index) > 0:
        df = df[['Sample Name', 'Index_ND', 'Species', 'Custom Reference']]
        df.to_csv(custom_refdata_csv_path, index=False)
        sys.stdout.write(f"custom_refdata.csv ({len(df.index)} lines written) \n")

    # check files exist
    assert os.path.exists(rna_counts_csv_path), 'ERROR: RNAcounts.tsv not created'
    assert os.path.exists(indexes_csv_path), 'ERROR: Indexes.csv not created'
    assert os.path.exists(sb_counts_csv_path), 'ERROR: SBcounts.tsv not created'
    assert os.path.exists(spatial_csv_path), 'ERROR: Spatial.tsv not created'
    assert os.path.exists(name_to_index_csv_path), 'ERROR: name_to_index.csv not created'
    
    return None


def main(bcl, src_path, data_path, default_bcl_path, sheet_id, reference_path):
    service = pars_google_sheet(src_path)
    rows, df = match_bcl_arg(bcl, service, sheet_id)
    get_info_and_create_files(rows, df, bcl, data_path, default_bcl_path, reference_path)

if __name__ == "__main__":
    main(bcl, src_path, data_path, default_bcl_path, sheet_id, reference_path)