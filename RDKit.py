# -*- coding: utf-8 -*-
"""
Calculate 200 2D RDKit descriptors
for either training or test sequence set

NOTE: Descriptors that are calculated depend on the
python version used. For replicability with this study,
Python 3.6.8 was used to calculated 200 descriptors for
each training or test sequence.

Author: Joshua D Richardson
"""

import numpy as np
import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

#############################################################################################################################
# BEGIN USER INPUT PARAMETERS

# Define these according to your local directory #
user = 'jdrichardso3'
home_dir = f'C:/Users/{user}/OneDrive - UW-Madison/Documents/Research/Papers/Paper Drafts/2_aB Peptide Predictions/Data/'

# Input parameters
data = 'test' # train test
test_bb = 'aaB' # IF data = 'test' --> aaB (Rounds 1-4) or aaBaaaB (Rounds 5-6)

# END USER INPUT PARAMETERS
##############################################################################################################################

# Training and Test Paths based on home_dir above
train_path = f'{home_dir}/Training'
test_path = f'{home_dir}/Test Design Space'

# Import Peptide Smiles csv file (use Make_Smiles.py first for 168,000 test sequence space)
if data == 'train':
    smi_df = pd.read_csv(f'{train_path}/SMILES.csv', index_col = 0, header = 0)
elif data == 'test':
    if test_bb == 'aaB':
        smi_df = pd.read_csv(f'{test_path}/SMILES_aaB.csv', index_col = 0, header = 0)
    elif test_bb == 'aaBaaaB':
        smi_df = pd.read_csv(f'{test_path}/SMILES_aaBaaaB.csv', index_col = 0, header = 0)

# Extract SMILES columns
idx = smi_df.index
SMILES = list(smi_df['SMILES'])

# Loop for calculating RDKit descriptors / initialize descriptor matrix
descriptors_list = [x[0] for x in Descriptors._descList]
descriptor_matrix = np.zeros((len(SMILES), len(descriptors_list)))

iteration = 1
for i in range(len(idx)):
    m1 = Chem.MolFromSmiles(SMILES[i])
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(descriptors_list)
    ds = calc.CalcDescriptors(m1)
    descriptor_matrix[i,:] = ds
    if iteration%1000 == 0:
        print(f'Done with sequence {iteration} of {len(SMILES)}')
    iteration+=1

# Save descriptor matrix to dataframe
desc_mat_df = pd.DataFrame(descriptor_matrix, columns = descriptors_list)
if data == 'train':
    desc_mat_df.to_csv(f'{train_path}/RDKit.csv', encoding='utf-8-sig')
elif data == 'test':
    if test_bb == 'aaB':
        desc_mat_df.to_csv(f'{train_path}/RDKit_aaB.csv', encoding='utf-8-sig')
    elif test_bb == 'aaBaaaB':
        desc_mat_df.to_csv(f'{train_path}/RDKit_aaBaaaB.csv', encoding='utf-8-sig')



    