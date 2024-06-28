import os
import gzip
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from argparse import ArgumentParser

# Calculate RMSD and CV
def rmsd(x):
    return np.sqrt(np.mean((x - x.mean())**2))

def cv(x):
    return np.std(x) / x.mean()

# Parse residue number
def parse_residue_number(residue):
    return int(''.join(filter(str.isdigit, residue)))

# Process a single file
def process_file(file_path):
    try: # Try to open and read the file in gzip format
        with gzip.open(file_path, 'rt') as f:
            df = pd.read_csv(f)
    except OSError: # If the file is not in gzip format
        with open(file_path, 'r') as f:
            df = pd.read_csv(f)
    
    tmp = df.groupby(['Residue1', 'Residue2', 'Range']).agg({
        'Surface': 'sum',
        'Volume' : 'sum',
        'AOWV'   : 'sum'
    }).reset_index()
    
    tmp['Residue1_last'] = tmp['Residue1'].str[-1] # Extract the last digit of Residue1
    tmp['Residue2_last'] = tmp['Residue2'].str[-1] # Extract the last digit of Residue2

    return tmp

# Create a parser
parser = ArgumentParser(description="Generate heatmap from data.")
# Add arguments
parser.add_argument('-i', nargs='?', default='out', help="Input directory.")
parser.add_argument('-col', nargs='*', default=['Surface_max', 'Surface_min', 'Surface_range',
           'Surface_mean', 'Surface_rmsd', 'Surface_cv',
           'Volume_max', 'Volume_min', 'Volume_range',
           'Volume_mean', 'Volume_rmsd', 'Volume_cv',
           'AOWV_max', 'AOWV_min', 'AOWV_range',
           'AOWV_mean', 'AOWV_rmsd', 'AOWV_cv'], help="Columns to be processed.")
# Parse arguments
args = parser.parse_args()

results = args.i
columns = args.col

try:
    i = 0
    for root, dirs, files in os.walk(results): # Iterate over all subfolders in the directory
        i += 1
        if i == 1: # Ensure only one iteration
            for dir in dirs:
                output_df = pd.DataFrame() # Create an empty DataFrame to store the results
                sub_dir = os.path.join(root, dir)
                for sub_root, sub_dirs, sub_files in os.walk(sub_dir):
                    for file in sub_files:
                        if file.startswith('_ALL'): # If the file name starts with 'ALL'
                            file_path = os.path.join(sub_root, file)
                            df = process_file(file_path)
                            output_df = pd.concat([output_df, df]) # Add the results to output_df
                if not output_df.empty: # If output_df is not empty
                    output_df['Pair'] = output_df['Residue1_last'] + '-' + output_df['Residue2_last']
                    output_df.drop(labels=['Residue1_last', 'Residue2_last'],axis=1,inplace=True)
                    
                    # Create a ID column, sort Residue1 and Residue2 by their numbers
                    output_df['ID'] = output_df.apply(lambda row: '-'.join(sorted([row['Residue1'], row['Residue2']], key=parse_residue_number)), axis=1)
                    # Group by ID, and calculate the Surface, Volume, and AOWV columns for each group
                    result = output_df.groupby('ID').filter(lambda x: len(x) >= 5).groupby('ID').agg({
                        'Pair': 'first',
                        'Surface' : ['max', 'min', lambda x: x.max() - x.min(), 'mean', rmsd, cv],
                        'Volume'  : ['max', 'min', lambda x: x.max() - x.min(), 'mean', rmsd, cv],
                        'AOWV': ['max', 'min', lambda x: x.max() - x.min(), 'mean', rmsd, cv],
                    }).reset_index()
                    # Rename the column names
                    result.columns = ['_'.join(col).strip() for col in result.columns.values]
                    result.rename(columns={'ID_': 'ID', 'Surface_<lambda_0>': 'Surface_range',
                                           'Volume_<lambda_0>': 'Volume_range',
                                           'AOWV_<lambda_0>': 'AOWV_range'}, inplace=True)
                    # Round the results to 3 decimal places
                    result = result.round(3)            
                    result.to_csv(f'{dir}.csv', index=False)  # Save the results to a csv file
                    
                    # Get all IDs and sort them by the residue number
                    ids = sorted(result['ID'].apply(lambda x: x.split('-')).explode().unique(), key=lambda x: int(x[1:-1]))

                    # Create an empty DataFrame to store the values for plotting
                    for column in columns:
                        col = pd.DataFrame(index=ids, columns=ids)
                        
                        try:
                            # Fill the DataFrame
                            for index, row in result.iterrows():
                                id1, id2 = row['ID'].split('-')
                                col.loc[id1, id2] = row[column]
                                col.loc[id2, id1] = row[column]

                            # Convert the data in the DataFrame to float
                            col = col.astype(float)
                            
                            # Create a heatmap
                            name = column.replace('_', ' ').title()
                            plt.figure(figsize=(30, 24), dpi=300)
                            sns.heatmap(col, annot=False, cmap='coolwarm', vmax=100)
                            plt.title(name)
                            plt.savefig(column+f'_{dir}.png')
                            plt.close()
                            
                        except:
                            continue
                                
except:
    print('The heatmap could not be successfully generated due to an error')
        
