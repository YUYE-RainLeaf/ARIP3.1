import math
from pathlib import Path

import numpy as np
import pandas as pd

from .utils import timer
from .typing import *


def range_kind(surface, volume, threshold):
    info_df = pd.DataFrame(list(surface.items()), columns=['Pair', 'Tuple'])
    info_df[['Surface', 'Distance', 'Type1', 'Type2']] = pd.DataFrame(info_df['Tuple'].tolist(), index=info_df.index)
    info_df[['Residue1', 'Atom1', 'Residue2', 'Atom2']] = info_df['Pair'].apply('+'.join).str.split('+', expand=True)
    
    # Add volume, atom type, interaction type
    info_df['Volume'] = np.nan
    info_df['Range']  = 0
    info_df['Kind']   = 'UNDEF'
    info_df['Dist1']  = info_df['Residue1'].apply(lambda x: int(x.split(',')[0].split(';')[0][1:]))
    info_df['Dist2']  = info_df['Residue2'].apply(lambda x: int(x.split(',')[0].split(';')[0][1:]))
    info_df['Dist']   = info_df['Dist1'] - info_df['Dist2']
    
    # Match volume
    if volume != {}:
        info_df['Volume'] = 0
        Volu = pd.DataFrame(list(volume.items()), columns=['Pair', 'Volume'])
        info_df = info_df.merge(Volu[['Pair', 'Volume']], on='Pair', how='left')
        info_df = info_df.rename(columns={'Volume_y': 'Volume'})
    
        # Some contacts with almost zero volume return null values during calculation, they are invalid data
        info_df = info_df.dropna(subset=['Volume'])
     
    # Match interaction type
        # Below are proteins
    info_df['Kind'] = np.where((info_df['Type1'] != 'D') & (info_df['Type1'] != 'R') & # Other interactions of amino acids
                               (info_df['Type1'] != 'X') & (info_df['Type2'] != 'D') &
                               (info_df['Type2'] != 'R') & (info_df['Type2'] != 'X'), 'OTHER', info_df['Kind'])
    
    info_df['Kind'] = np.where((info_df['Type1'] == 'I') & # Hydrogen bond
                               (info_df['Type2'] == 'I') &
                               (1.5 <= info_df['Distance']) &
                               (info_df['Distance'] <= 3.5), 'HB', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'I') &
                               (info_df['Type2'] == 'II') &
                               (1.5 <= info_df['Distance']) &
                               (info_df['Distance'] <= 3.5), 'HB', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'I') &
                               (info_df['Type2'] == 'III') &
                               (1.5 <= info_df['Distance']) &
                               (info_df['Distance'] <= 3.5), 'HB', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'II') &
                               (info_df['Type2'] == 'I') &
                               (1.5 <= info_df['Distance']) &
                               (info_df['Distance'] <= 3.5), 'HB', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'II') &
                               (info_df['Type2'] == 'III') &
                               (1.5 <= info_df['Distance']) &
                               (info_df['Distance'] <= 3.5), 'HB', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'III') &
                               (info_df['Type2'] == 'I') &
                               (1.5 <= info_df['Distance']) &
                               (info_df['Distance'] <= 3.5), 'HB', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'III') &
                               (info_df['Type2'] == 'II') &
                               (1.5 <= info_df['Distance']) &
                               (info_df['Distance'] <= 3.5), 'HB', info_df['Kind'])
    
    info_df['Kind'] = np.where((info_df['Type1'] == 'V') & # Aromatic interaction
                               (info_df['Type2'] == 'V'), 'AROM', info_df['Kind'])
    
    info_df['Kind'] = np.where((info_df['Type1'] == 'IV') & # Hydrophobic interaction
                               (info_df['Type2'] == 'IV'), 'PHOB', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'IV') &
                               (info_df['Type2'] == 'V'),  'PHOB', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'V') &
                               (info_df['Type2'] == 'IV'), 'PHOB', info_df['Kind'])
    
    info_df['Kind'] = np.where((info_df['Type1'] == 'I') & # Destabilizing interaction
                               (info_df['Type2'] == 'IV'), 'DC', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'II') &
                               (info_df['Type2'] == 'II'), 'DC', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'II') &
                               (info_df['Type2'] == 'IV'), 'DC', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'II') &
                               (info_df['Type2'] == 'VIII'), 'DC', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'III') &
                               (info_df['Type2'] == 'III'), 'DC', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'III') &
                               (info_df['Type2'] == 'IV'), 'DC', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'III') &
                               (info_df['Type2'] == 'VII'), 'DC', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'IV') &
                               (info_df['Type2'] == 'I'), 'DC', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'IV') &
                               (info_df['Type2'] == 'II'), 'DC', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'IV') &
                               (info_df['Type2'] == 'III'), 'DC', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'VII') &
                               (info_df['Type2'] == 'III'), 'DC', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'VII') &
                               (info_df['Type2'] == 'VII'), 'DC', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'VIII') &
                               (info_df['Type2'] == 'II'), 'DC', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'VIII') &
                               (info_df['Type2'] == 'VIII'), 'DC', info_df['Kind'])
    
    info_df['Kind'] = np.where((info_df['Atom1'] == 'C') & # Covalent peptide bond, and the next residue
                               (info_df['Atom2'] == 'N') &
                               (info_df['Type1'] != 'X') & (info_df['Type2'] != 'X') &
                               ((info_df['Residue1'].apply(lambda x: x[0]) ==
                                 info_df['Residue2'].apply(lambda x: x[0]))) &
                               (info_df['Dist'] == -1), 'PB', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Atom1'] == 'N') & # Covalent peptide bond, and the previous residue
                               (info_df['Atom2'] == 'C') &
                               (info_df['Type1'] != 'X') & (info_df['Type2'] != 'X') &
                               ((info_df['Residue1'].apply(lambda x: x[0]) ==
                                 info_df['Residue2'].apply(lambda x: x[0]))) &
                               (info_df['Dist'] == 1), 'PB', info_df['Kind'])
    
    info_df['Kind'] = np.where((info_df['Atom1'] == 'SG') & # Covalent disulfide bond
                               (info_df['Atom2'] == 'SG') &
                               (1.95 <= info_df['Distance']) &
                               (info_df['Distance'] <= 2.1), 'SS', info_df['Kind'])
    
        # Below are nucleic acids
    info_df['Kind'] = np.where((info_df['Type1'] == 'D') & # Interactions between deoxyribonucleotides
                               (info_df['Type2'] == 'D'), 'DD', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'R') & # Interactions between ribonucleotides
                               (info_df['Type2'] == 'R'), 'RR', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'D') & # Interactions between ribonucleotides and deoxyribonucleotides
                               (info_df['Type2'] == 'R'), 'DR', info_df['Kind'])    
    info_df['Kind'] = np.where((info_df['Type1'] == 'R') & # Interactions between ribonucleotides and deoxyribonucleotides
                               (info_df['Type2'] == 'D'), 'DR', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] != 'D') & (info_df['Type1'] != 'R') & # Interactions between amino acids and deoxyribonucleotides
                               (info_df['Type1'] != 'X') & (info_df['Type2'] == 'D'), 'AD', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'D') & (info_df['Type2'] != 'R') & # Interactions between amino acids and deoxyribonucleotides
                               (info_df['Type2'] != 'X') & (info_df['Type2'] != 'D'), 'AD', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] != 'D') & (info_df['Type1'] != 'R') & # Interactions between amino acids and ribonucleotides
                               (info_df['Type1'] != 'X') & (info_df['Type2'] == 'R'), 'AR', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Type1'] == 'R') & (info_df['Type2'] != 'R') & # Interactions between amino acids and ribonucleotides
                               (info_df['Type2'] != 'X') & (info_df['Type2'] != 'D'), 'AR', info_df['Kind'])

    info_df['Kind'] = np.where((info_df['Atom1'] == "O3'") & # Covalent phosphodiester bond, and the next residue
                               (info_df['Atom2'] == 'P') & (info_df['Type2'] != 'X') &
                               ((info_df['Residue1'].apply(lambda x: x[0]) ==
                                 info_df['Residue2'].apply(lambda x: x[0]))) &
                               (info_df['Dist'] == -1), 'PD', info_df['Kind'])
    info_df['Kind'] = np.where((info_df['Atom1'] == 'P') & # Covalent phosphodiester bond, and the previous residue
                               (info_df['Atom2'] == "O3'") & (info_df['Type1'] != 'X') &
                               ((info_df['Residue1'].apply(lambda x: x[0]) ==
                                 info_df['Residue2'].apply(lambda x: x[0]))) &
                               (info_df['Dist'] == 1), 'PD', info_df['Kind'])
    
    # Match interaction distance range
    info_df['Range'] = np.where((abs(info_df['Dist']) <= 2) & # Short range
                                ((info_df['Residue1'].apply(lambda x: x[0]) ==
                                  info_df['Residue2'].apply(lambda x: x[0]))), 'S', info_df['Range'])
    info_df['Range'] = np.where((abs(info_df['Dist']) > 2) & # Medium range
                                (abs(info_df['Dist']) <= 4) &
                                ((info_df['Residue1'].apply(lambda x: x[0]) ==
                                  info_df['Residue2'].apply(lambda x: x[0]))), 'M', info_df['Range'])
    info_df['Range'] = np.where((abs(info_df['Dist']) > 4) & # Long range
                                (abs(info_df['Dist']) <= 10) &
                                ((info_df['Residue1'].apply(lambda x: x[0]) ==
                                  info_df['Residue2'].apply(lambda x: x[0]))), 'L1', info_df['Range'])
    info_df['Range'] = np.where((abs(info_df['Dist']) > 10) &
                                (abs(info_df['Dist']) <= 20) &
                                ((info_df['Residue1'].apply(lambda x: x[0]) ==
                                  info_df['Residue2'].apply(lambda x: x[0]))), 'L2', info_df['Range'])
    info_df['Range'] = np.where((abs(info_df['Dist']) > 20) &
                                (abs(info_df['Dist']) <= 30) &
                                ((info_df['Residue1'].apply(lambda x: x[0]) ==
                                  info_df['Residue2'].apply(lambda x: x[0]))), 'L3', info_df['Range'])
    info_df['Range'] = np.where((abs(info_df['Dist']) > 30) &
                                (abs(info_df['Dist']) <= 40) &
                                ((info_df['Residue1'].apply(lambda x: x[0]) ==
                                  info_df['Residue2'].apply(lambda x: x[0]))), 'L4', info_df['Range'])
    info_df['Range'] = np.where((abs(info_df['Dist']) > 40) &
                                (abs(info_df['Dist']) <= 50) &
                                ((info_df['Residue1'].apply(lambda x: x[0]) ==
                                  info_df['Residue2'].apply(lambda x: x[0]))), 'L5', info_df['Range'])
    info_df['Range'] = np.where((abs(info_df['Dist']) > 50) &
                                ((info_df['Residue1'].apply(lambda x: x[0]) ==
                                  info_df['Residue2'].apply(lambda x: x[0]))), 'L6', info_df['Range'])
    # Interchain
    info_df['Range'] = np.where((info_df['Residue1'].apply(lambda x: x[0]) !=
                                 info_df['Residue2'].apply(lambda x: x[0])), 'I', info_df['Range'])
    
    # Volume might be a float, but if the volume based on electron density is calculated, it will be a list
    is_list = info_df['Volume'].apply(lambda x: isinstance(x, list))
    if is_list.any() == True:
        tmp_df = info_df.loc[is_list, 'Volume'].apply(pd.Series)
        tmp_df.columns = ['Volume', 'AOWV']
        info_df = pd.concat([info_df.drop(['Volume'], axis=1), tmp_df], axis=1)

        all_df = info_df[['Residue1', 'Atom1', 'Type1', 'Residue2', 'Atom2', 'Type2',
                          'Distance', 'Surface', 'Volume', 'AOWV', 'Range', 'Kind']]
    
    else:
        all_df = info_df[['Residue1', 'Atom1', 'Type1', 'Residue2', 'Atom2', 'Type2',
                          'Distance', 'Surface', 'Volume', 'Range', 'Kind']]
    
    # If screen by threshold
    if threshold is not None:
        if volume != {}:
            all_df = all_df[(all_df['Surface'] >= threshold[0]) & (all_df['Volume'] >= threshold[1])]
        else:
            all_df = all_df[(all_df['Surface'] >= threshold[0])]
    
    # Sort by chain and sequence
    all_df_copy = all_df.copy()
    all_df_copy['Chain1']    = all_df_copy['Residue1'].str[0]
    all_df_copy['Sequence1'] = all_df_copy['Residue1'].apply(lambda x: int(x.split(',')[0].split(';')[0][1:]))    
    all_df_copy['Chain2']    = all_df_copy['Residue2'].str[0]
    all_df_copy['Sequence2'] = all_df_copy['Residue2'].apply(lambda x: int(x.split(',')[0].split(';')[0][1:]))
    all_df_copy = all_df_copy.sort_values(['Chain1', 'Sequence1', 'Chain2', 'Sequence2'])

    all_df = all_df_copy.drop(columns=['Chain1', 'Sequence1', 'Chain2', 'Sequence2'])

    return all_df


def summary(sasa, dihedral_angle, Residues):
    # Summary of each residue
    EMPTY = np.zeros([len(Residues)])
    sum_df = pd.DataFrame()
    sum_df['Residue']    = list(Residues.keys())
    sum_df['Phi']        = EMPTY.copy()
    sum_df['Psi']        = EMPTY.copy()
    sum_df['SASA']       = np.nan # SASA, add first, if not available, delete later
    sum_df['Cova_Surf']  = EMPTY.copy()
    sum_df['Cova_Volu']  = np.nan
    sum_df['Cova_AOWV']  = np.nan # Add first, if not available, delete later
    sum_df['NC_Surf']    = EMPTY.copy()
    sum_df['NC_Volu']    = np.nan
    sum_df['NC_AOWV']    = np.nan # Add first, if not available, delete later
    sum_df['UNDEF_Surf'] = np.nan 
    sum_df['UNDEF_Volu'] = np.nan 
    sum_df['UNDEF_AOWV'] = np.nan # Unknown type of interaction, add first, if not available, delete later
    
    for i in range(len(sum_df)):
        a_res = sum_df.loc[i]

        # Match SASA
        if sasa:
            if a_res['Residue'] in sasa:
                sa = sasa[a_res['Residue']]
                sum_df.loc[i, 'SASA'] = round(math.degrees(sa), 3)

        # Match dihedral angle
        if dihedral_angle:
            if a_res['Residue'] in dihedral_angle:
                phi, psi = dihedral_angle[a_res['Residue']]
            else:
                phi, psi = None, None
        else:
            phi, psi = None, None
        
        # Convert to degrees or set to Nan
        sum_df.loc[i, 'Psi'] = np.nan if psi is None else round(math.degrees(psi), 3)
        sum_df.loc[i, 'Phi'] = np.nan if phi is None else round(math.degrees(phi), 3)
        
        # Non-covalent and covalent interaction surface and volume
        interaction = Residues[a_res['Residue']]
        NC   = interaction[~interaction['Kind'].isin(['SS', 'PB', 'PD', 'UNDEF'])]
        Cova = interaction[ interaction['Kind'].isin(['SS', 'PB', 'PD'])]
        UD   = interaction[ interaction['Kind'].isin(['UNDEF'])]
        
        # Surf
        sum_df.loc[i, 'Cova_Surf'] = np.nansum(Cova['Surface'])
        sum_df.loc[i, 'NC_Surf']   = np.nansum(NC['Surface'])
        if 'Volume' in interaction.columns:
            sum_df.loc[i, 'Cova_Volu'] = np.nansum(Cova['Volume'])
            sum_df.loc[i, 'NC_Volu']   = np.nansum(NC['Volume'])
        
        if 'AOWV' in interaction.columns:
            sum_df.loc[i, 'Cova_AOWV'] = np.nansum(Cova['AOWV'])
            sum_df.loc[i, 'NC_AOWV']   = np.nansum(NC['AOWV'])
        
        if not UD.empty:
            sum_df.loc[i, 'UNDEF_Surf'] = np.nansum(UD['Surface'])
            if 'Volume' in interaction.columns:
                sum_df.loc[i, 'UNDEF_Volu'] = np.nansum(UD['Volume'])
            if 'AOWV' in interaction.columns:
                sum_df.loc[i, 'UNDEF_AOWV']  = np.nansum(UD['AOWV'])
           
    # Round the columns 'Cova_Surf', 'Cova_Volu', 'Cova_AOWV', 'NC_Surf', 'NC_Volu', 'NC_AOWV',
    # 'UNDEF_Surf', 'UNDEF_Volu', 'UNDEF_AOWV' to three decimal places
    sum_df[['Cova_Surf',  'Cova_Volu',  'Cova_AOWV',
            'NC_Surf',    'NC_Volu',    'NC_AOWV',
            'UNDEF_Surf', 'UNDEF_Volu', 'UNDEF_AOWV']] = sum_df[['Cova_Surf',  'Cova_Volu',  'Cova_AOWV',
                                                                'NC_Surf',    'NC_Volu',    'NC_AOWV',
                                                                'UNDEF_Surf', 'UNDEF_Volu', 'UNDEF_AOWV']].round(3)
    
    # Sort by chain and sequence
    sum_df['Chain']    = sum_df['Residue'].str[0]
    sum_df['Sequence'] = sum_df['Residue'].apply(lambda x: int(x.split(',')[0].split(';')[0][1:]))
    sorted_df = sum_df.sort_values(['Chain', 'Sequence'])

    # If all are nucleic acids, there are no dihedral angle parameters
    if set(sorted_df['Phi']) == {''} and set(sorted_df['Phi']) == {''}: # If both are empty values, then they are two 'True'
        res_df = sorted_df.drop(columns=['Chain', 'Sequence', 'Phi', 'Psi']).dropna(axis=1, how='all')
    else:
        res_df = sorted_df.drop(columns=['Chain', 'Sequence']).dropna(axis=1, how='all')

    res_df['Residue'] = res_df['Residue'].str.replace(',', '')
    res_df['Residue'] = res_df['Residue'].str.replace(';', '+')

    # Drop columns that are all zeros
    for col in res_df.columns:
        if (res_df[col] == 0).all():
            res_df.drop(col, axis=1, inplace=True)

    return res_df


def residue_pairs(Residues):
    final_rows = []

    for residue1, df in Residues.items():
        # Check for existence
        has_volu = 'Volume' in df.columns
        has_aowv = 'AOWV'   in df.columns
        
        # Group by Residue2 and Kind, and sum Surface and Volume
        grouped = df.groupby(['Residue2', 'Kind']).agg({'Range': 'first', 'Surface': 'sum',
                                                        **({'Volume': 'sum'} if has_volu else {}),
                                                        **({'AOWV': 'sum'} if has_aowv else {})}).reset_index()
       
        for _, row in grouped.iterrows():
            kind = row['Kind']
            prefix = ''
            
            # Determine the prefix
            if kind in ['SS', 'PB', 'PD']:
                prefix = 'Cova_'
            elif kind == 'UNDEF':
                prefix = 'UNDEF_'
            else:
                prefix = 'NC_'
            
            current_row = {
                'Residue1': residue1,
                'Residue2': row['Residue2'],
                'Range'   : row['Range'],
                'Cova_Surf' : 0     , 'Cova_Volu' : 0     , 'Cova_AOWV' : np.nan,
                'NC_Surf'   : 0     , 'NC_Volu'   : 0     , 'NC_AOWV'   : np.nan,
                'UNDEF_Surf': np.nan, 'UNDEF_Volu': np.nan, 'UNDEF_AOWV': np.nan,
            }
            
            current_row[prefix + 'Surf'] = row['Surface']
            if has_volu:
                current_row[prefix + 'Volu'] = row['Volume']
            if has_aowv:
                current_row[prefix + 'AOWV'] = row['AOWV']
        
            final_rows.append(current_row)
        
    final_df = pd.DataFrame(final_rows)
    aggregated_df = final_df.groupby(['Residue1', 'Residue2']).agg({'Range': 'first',
        'Cova_Surf' : 'sum', 'Cova_Volu' : 'sum', 'Cova_AOWV' : 'sum',
        'NC_Surf'   : 'sum', 'NC_Volu'   : 'sum', 'NC_AOWV'   : 'sum',
        'UNDEF_Surf': 'sum', 'UNDEF_Volu': 'sum', 'UNDEF_AOWV': 'sum',
    }).reset_index()

    aggregated_df.dropna(axis=1, how='all', inplace=True)
    # Sort
    def custom_sort_key(x):
        # Split the string to obtain the chain id and number
        parts = x.split(',') if ',' in x else x.split(';')
        chain_part = parts[0]
        chain = chain_part[0]
        number = int(chain_part[1:])
        return (chain, number)

    aggregated_df.sort_values(by=['Residue1', 'Residue2'], key=lambda x: x.map(custom_sort_key), inplace=True)
    aggregated_df['Residue1'] = aggregated_df['Residue1'].str.replace(',', '').str.replace(';', '+')
    aggregated_df['Residue2'] = aggregated_df['Residue2'].str.replace(',', '').str.replace(';', '+')    

    # Drop columns that are all zeros
    for col in aggregated_df.columns:
        if (aggregated_df[col] == 0).all():
            aggregated_df.drop(col, axis=1, inplace=True)
        
    return aggregated_df


def pdb_csv_sort(idx:int, name:str, dihedral_angle:Angles, sasa:SASA, surface:Surfaces, volume:Volumes, out_path:Path, threshold:List[float], compress:bool, each:bool, disable_print=False):
    @timer(disable_print=disable_print)
    def make_result():
        out_dp = out_path / (name+'_z') if compress else out_path / name
        out_dp.mkdir(parents=True, exist_ok=True)
        
        all_df = range_kind(surface, volume, threshold)
        all_df = all_df.dropna(axis=1, how='all')
        
        all_tmp = all_df.copy()
        all_tmp['Residue1'] = all_tmp['Residue1'].str.replace(',', '')
        all_tmp['Residue1'] = all_tmp['Residue1'].str.replace(';', '+')
        all_tmp['Residue2'] = all_tmp['Residue2'].str.replace(',', '')
        all_tmp['Residue2'] = all_tmp['Residue2'].str.replace(';', '+')

        Residues = dict(list(all_df.groupby('Residue1')))            
        res_df   = summary(sasa, dihedral_angle, Residues)
        pair_df  = residue_pairs(Residues)            
        
        # If there is only one MODEL
        if idx == -1:
            if compress:
                all_tmp.to_csv(out_dp / f'_ALL_{name}.gz', compression='gzip', index=0)
                res_df .to_csv(out_dp / f'_SUM_{name}.gz', compression='gzip', index=0)
                pair_df.to_csv(out_dp / f'_RES_{name}.gz', compression='gzip', index=0)
                if each:
                    Results  = dict(list(all_tmp.groupby('Residue1')))
                    for res in Results:
                        Results[res].to_csv(out_dp / f'{res}.gz', compression='gzip', index=0)

            else:
                all_tmp.to_csv(out_dp / f'_ALL_{name}.csv', index=0)    # Summary of contacts
                res_df .to_csv(out_dp / f'_SUM_{name}.csv', index=0)    # Dihedral angles, etc.
                pair_df.to_csv(out_dp / f'_RES_{name}.csv', index=0)    # Residue interactions
                if each:
                    Results  = dict(list(all_tmp.groupby('Residue1')))
                    for res in Results:
                        Results[res].to_csv(out_dp / f'{res}.csv', index=0)
                        
        # If there are multiple MODELs
        elif idx >= 0:
            model_dp = out_dp / f'{name}_{idx}'
            model_dp.mkdir(parents=True, exist_ok=True)
            
            if compress:
                all_tmp.to_csv(model_dp / f'_ALL_{name}.gz', compression='gzip', index=0)
                res_df .to_csv(model_dp / f'_SUM_{name}.gz', compression='gzip', index=0)
                pair_df.to_csv(model_dp / f'_RES_{name}.gz', compression='gzip', index=0)
                if each:
                    Results  = dict(list(all_tmp.groupby('Residue1')))
                    for res in Results:
                        Results[res].to_csv(model_dp / f'{res}.gz', compression='gzip', index=0)

            else:
                all_tmp.to_csv(model_dp / f'_ALL_{name}.csv', index=0)
                res_df .to_csv(model_dp / f'_SUM_{name}.csv', index=0)
                pair_df.to_csv(model_dp / f'_RES_{name}.csv', index=0)
                if each:
                    Results  = dict(list(all_tmp.groupby('Residue1')))
                    for res in Results:
                        Results[res].to_csv(model_dp / f'{res}.csv', index=0)

    make_result()
    
    