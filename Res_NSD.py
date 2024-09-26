# -*- coding: utf-8 -*-
"""
Calculate the NSD per round for newly
introduced residues (Nle, Q, Nva, βY, W)
in prediction rounds 1-4

Author: Joshua D Richardson
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

#############################################################################################################################
# BEGIN USER INPUT PARAMETERS

# Define these according to your local directory #
user = 'jdrichardso3'
home_dir = f'C:/Users/{user}/OneDrive - UW-Madison/Documents/Research/Papers/Paper Drafts/2_aB Peptide Predictions/Data/'

# Inputs
label = 'HC10' # HC10 MIC
residue = 'Gln' # Nle Gln (Q) Nva βTyr (βY) Trp (W)  <-- use 3 letter code for residue specification

# END USER INPUT PARAMETERS
##############################################################################################################################

# Training and Test Paths based on home_dir above
train_path = f'{home_dir}/Training'
test_path = f'{home_dir}/Test Design Space'

# Graphing parameters
plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=16)          # controls default text sizes
plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=16)    # legend fontsize
plt.rc('figure', titlesize=16)  # fontsize of the figure title

# Set up NSD dicts
nsd_dict = dict.fromkeys(['1', '2', '3', '4*', '4'])
seq_dict = dict.fromkeys(['1', '2', '3', '4*', '4'])

# First 4 rounds used in study
for idx, round_num in enumerate(['1st', '2nd', '3rd', '4th']):
    sd_res = []
    sd_seq = []
    test_df = pd.read_csv(f'{test_path}/Predictions/Reduced Test/Total/Results_{round_num}_Round.csv', header=0).set_index('Unnamed: 0', drop=True)
    for i in range(len(test_df)):
        seq_arr = test_df.iloc[i]['Sequence'].split('-')
        if residue in seq_arr:
            sd_res = np.append(sd_res,test_df.iloc[i][f'Log2{label} s.d.'])
            sd_seq = np.append(sd_seq,test_df.iloc[i]['Sequence'])
    print(f'Done with Round {idx+1}')
        
    norm_sd_res = [(i - min(test_df[f'Log2{label} s.d.']))/(max(test_df[f'Log2{label} s.d.'])-min(test_df[f'Log2{label} s.d.'])) for i in sd_res]
    nsd_dict[str(idx+1)] = norm_sd_res
    seq_dict[str(idx+1)] = sd_seq
    
# Additional hypothetical 4th Round (4*) where descriptors are kept the same as Rounds 1-3
sd_res = []
sd_seq = []

test_df = pd.read_csv(f'{test_path}/Predictions/Reduced Test/Total/Results_4th_Round_noupdate.csv', header=0).set_index('Unnamed: 0', drop=True)
for i in range(len(test_df)):
    seq_arr = test_df.iloc[i]['Sequence'].split('-')
    if residue in seq_arr:
        sd_res = np.append(sd_res,test_df.iloc[i][f'Log2{label} s.d.'])
        sd_seq = np.append(sd_seq,test_df.iloc[i]['Sequence'])
        
print('Done with Round 4*')

norm_sd_res = [(i - min(test_df[f'Log2{label} s.d.']))/(max(test_df[f'Log2{label} s.d.'])-min(test_df[f'Log2{label} s.d.'])) for i in sd_res]
nsd_dict['4*'] = norm_sd_res
seq_dict['4*'] = sd_seq

# Convert raw NSD vs Round data to DataFrame and save as .csv
nsd_df = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in nsd_dict.items()]))
nsd_df.to_csv(f'{test_path}/NSD Per Round/{label}_{residue}.csv')

# Plot bar and whisker plot and save
fig, ax = plt.subplots()

ax.set_box_aspect(1)
ax.set_xlabel('Prediction Round')
ax.set_ylabel('NSD')

data = [nsd_dict['1'], nsd_dict['2'], nsd_dict['3'], nsd_dict['4*'], nsd_dict['4']]
meanpointprops = dict(marker='s', markeredgecolor='black',
                      markerfacecolor='black')
bp = ax.boxplot(data, whis=(10,90), showfliers='no', sym="", showmeans=True, meanprops=meanpointprops)

plt.xticks([1, 2, 3, 4, 5], ['1', '2', '3', '4*', '4'])
plt.yticks([0,0.5,1])
plt.ylim([0,1])

fig.savefig(f'{test_path}/Figures/NSD Per Round/{label}_{residue}.svg', format='svg', bbox_inches='tight')
plt.show()

