# -*- coding: utf-8 -*-
"""
Script to find how many residues different
test sequences are from any training sequence

('Dif'' Column in Figure 4a)

Author: Joshua D Richardson
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

#############################################################################################################################
# BEGIN USER INPUT PARAMETERS

# Define these according to your local directory #
user = 'jdrichardso3'
home_dir = f'C:/Users/{user}/OneDrive - UW-Madison/Documents/Research/Papers/Paper Drafts/2_aB Peptide Predictions/Data/'

# Input parameters
round_num = '1st' # 1st 2nd 3rd 4th 5th 6th

# END USER INPUT PARAMETERS
##############################################################################################################################

# Training and Test Paths based on home_dir above
train_path = f'{home_dir}/Training'
test_path = f'{home_dir}/Test Design Space'

# How many sequences available per round
if round_num == '1st':
    row_keep = 147
elif round_num == '2nd':
    row_keep = 150
elif round_num == '3rd':
    row_keep = 154 
elif round_num == '4th':
    row_keep = 158
elif round_num == '5th':
    row_keep = 162
elif round_num == '6th':
    row_keep = 166

# Import sequences training peptides in round defined above
df = pd.read_excel(f'{train_path}/Sequence_Activity.xlsx', header=0)
df = df[:row_keep]
pep_num = df['Index'].dropna()
seq = df['Sequence']

# Backbone array
bb = df['Backbone']
    
# Find indices where bb = aaB (Rounds 1-4) or aaBaaaB (Rounds 5-6)
bb_idx = np.array([])
for idx, i in enumerate(bb):
    if round_num=='1st' or round_num=='2nd' or round_num=='3rd' or round_num=='4th':
        if i == 'ααβ':
            bb_idx = np.append(bb_idx, idx)
    elif round_num=='5th' or round_num=='6th':
        if i == 'ααβαααβ':
            bb_idx = np.append(bb_idx, idx)
bb_ref = bb.iloc[bb_idx]

# Relate bb_ref to seq of same aaB or aaBaaaB backbone
seq_ref = seq.iloc[bb_idx]
seq_ref_mat = seq_ref.str.split('-',expand=True)

# Initialize num_res_dif array
num_res_dif = np.array([])

# Import sequences of 168,000 new peptides
if round_num=='1st' or round_num=='2nd' or round_num=='3rd' or round_num=='4th':
    new_df = pd.read_csv(f'{test_path}/SMILES_aaB.csv', header=0, index_col=0)  # aaB backbone for Rounds 1-4
else:
    new_df = pd.read_csv(f'{test_path}/SMILES_aaBaaaB.csv', header=0, index_col=0)  # aaBaaaB backbone for Rounds 5-6
new_seq = new_df['seq']
new_seq_mat = new_seq.str.split('-',expand=True)

# Calculate peptide test seq (new_seq_mat) that are 1,2,3,4 seq different from any training (seq_ref_mat)
if os.path.exists(f'{test_path}/Res Dif/{round_num}_Round.csv') == False:
    
    # Loop to calculate peptide seq that are 1,2,3,4 seq different from each other
    within_one_idx = []
    for i in range(len(new_seq_mat.index)):
        sim_aa = np.array([])
        for j in range(len(seq_ref_mat.index)):
            sim_aa = np.append(sim_aa, np.count_nonzero(new_seq_mat.iloc[i]==seq_ref_mat.iloc[j]))
            
        # Conditionals based on residue differences
        if max(sim_aa) == 14:
            num_res_dif = np.append(num_res_dif,0)
        elif max(sim_aa) == 13:
            num_res_dif = np.append(num_res_dif,1)
        elif max(sim_aa) == 12:
            num_res_dif = np.append(num_res_dif,2)
        elif max(sim_aa) == 11:
            num_res_dif = np.append(num_res_dif,3)
        elif max(sim_aa) == 10:
            num_res_dif = np.append(num_res_dif,4)
            
        if i%1000 == 0:
            print(i)
    
    num_res_dif_df = pd.DataFrame()
    num_res_dif_df['Dif'] = num_res_dif
    num_res_dif_df.to_csv(f'{test_path}/Res Dif/{round_num}_Round.csv', encoding='utf-8-sig')
    
else:
    num_res_dif_df = pd.read_csv(f'{test_path}/Res Dif/{round_num}_Round.csv', index_col = 0, header = 0)

# Histogram to visualize # residues different from training
fig, ax = plt.subplots()
ax.set_box_aspect(1)
#ax.set_title('# Residues Different from Training')
ax.set_xlabel('# Residues Different')
ax.set_ylabel('Count')
ax.set_ylim([0,168000])
bins = np.arange(0, 4 + 1.5) - 0.5
counts, edges, bars = ax.hist(num_res_dif_df['Test Index'], bins = bins, rwidth=0.9)
ax.set_xticks(bins + 0.5)
plt.bar_label(bars)
plt.show()
