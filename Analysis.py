# -*- coding: utf-8 -*-
"""
Various process scripts useful for peptide
test sequence selection per round and
visualization of train/test distributions

Author: Joshua D Richardson
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#############################################################################################################################
# BEGIN USER INPUT PARAMETERS

# Define these according to your local directory #
user = 'jdrichardso3'
home_dir = f'C:/Users/{user}/OneDrive - UW-Madison/Documents/Research/Papers/Paper Drafts/2_aB Peptide Predictions/Data/'

# Inputs
round_num = '6th' # 1st 2nd 3rd 4th 5th 6th
reduced_test = 'yes' # Whether to look at all 168,000 test sequences ('no') or just reduced test space in Table S5 ('yes')

# What to plot
bb_train = 'yes' # Figure 1c-d & Figure S1
bb_type = 'ααβ' # If bb_train = yes (ααβ, ααβαααβ, αααβ, αβαβααβ)

design_space_pred = 'yes' # Figure S15 Log2(HC10) vs Log2(MIC) plots
design_space_sd = 'yes' # Figure S15 st. dev. distributions

# END USER INPUT PARAMETERS
##############################################################################################################################

# Training and Test Paths based on home_dir above
train_path = f'{home_dir}/Training'
test_path = f'{home_dir}/Test Design Space'

# Graphing Parameters
plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=16)          # controls default text sizes
plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=16)    # legend fontsize
plt.rc('figure', titlesize=16)  # fontsize of the figure title

# Number of sequences available per round
if round_num == '1st':
    row_keep = 147
if round_num == '2nd':
    row_keep = 150
elif round_num == '3rd':
    row_keep = 154 
elif round_num == '4th':
    row_keep = 158
elif round_num == '5th':
    row_keep = 162
elif round_num == '6th':
    row_keep = 166
elif round_num == '7th':
    row_keep = 170

# SI bands (Grey lines in Log2HC10 vs. Log2MIC plots)
log2MIC = [0.5,1,2,3,4,5,6,7,8,9,10]
log2HC10_SI10 = [np.log2(10*2**a) for a in log2MIC] # SI=10
log2HC10_SI20 = [np.log2(20*2**a) for a in log2MIC] # SI=20
log2HC10_SI30 = [np.log2(30*2**a) for a in log2MIC] # SI=30
log2HC10_SI40 = [np.log2(40*2**a) for a in log2MIC] # SI=40
log2HC10_SI50 = [np.log2(50*2**a) for a in log2MIC] # SI=50
 
##################
#### Training ####
##################

df = pd.read_excel(f'{train_path}/Sequence_Activity.xlsx', header=0)
df = df.head(row_keep)
pep_num = df['Index']
seq = df['Sequence']

y_train_HC10_raw = df['HC10']
y_train_HC10 = np.log2(y_train_HC10_raw)
y_train_MIC_raw = df['MIC']
y_train_MIC = np.log2(y_train_MIC_raw)

train_df = pd.DataFrame()
train_df['Sequence'] = list(seq)
train_df['Log2HC10'] = y_train_HC10
train_df['Log2MIC'] = list(y_train_MIC)
train_df['SI'] = 2**train_df['Log2HC10'] / 2**train_df['Log2MIC']

##################
####   Test   ####
##################

if reduced_test == 'yes':
    y_test_HC10 = pd.read_csv(f'{test_path}/Predictions/Reduced Test/HC10/HC10_Results_{round_num}_Round.csv', header=0).set_index('Unnamed: 0', drop=True)
    y_test_MIC = pd.read_csv(f'{test_path}/Predictions/Reduced Test/MIC/MIC_Results_{round_num}_Round.csv', header=0).set_index('Unnamed: 0', drop=True)
    intersect = set(y_test_HC10.index).intersection(set(y_test_MIC.index)) # Find index overlap between HC10/MIC predictions ('Total' column in Table S5)
    intersect = sorted(list(intersect))
    test_df = pd.DataFrame(index = intersect)
    
elif reduced_test == 'no':
    y_test_HC10 = pd.read_csv(f'{test_path}/Predictions/All 168000/HC10/HC10_Results_{round_num}_Round.csv', header=0).set_index('Unnamed: 0', drop=True)
    y_test_MIC = pd.read_csv(f'{test_path}/Predictions/All 168000/MIC/MIC_Results_{round_num}_Round.csv', header=0).set_index('Unnamed: 0', drop=True)
    test_df = pd.DataFrame()

# Import test sequences to add as column in test_df
if round_num == '1st' or round_num == '2nd' or round_num == '3rd' or round_num == '4th':
    new_df = pd.read_csv(f'{test_path}/SMILES_aaB.csv', header =0 , index_col=0)
else:
    new_df = pd.read_csv(f'{test_path}/SMILES_aaBaaaB.csv', header =0 , index_col=0)
new_seq = new_df['Sequence']
test_df['Sequence'] = new_seq

# Add # residues different from training as column if available (from Res_Dif.py)
if os.path.exists(f'{test_path}/Res Dif/{round_num}_Round.csv') == True:
    num_res_dif_df = pd.read_csv(f'{test_path}/Res Dif/{round_num}_Round.csv', index_col = 0, header = 0)
    num_res_dif = num_res_dif_df['Dif']
    test_df['Dif'] = num_res_dif
    
# Add HC10 and MIC predictions
test_df['Log2HC10'] = y_test_HC10['HC10']
test_df['Log2HC10 s.d.'] = y_test_HC10['Stdev']
test_df['Log2MIC'] = y_test_MIC['MIC']
test_df['Log2MIC s.d.'] = y_test_MIC['Stdev']

# Add SI (Calculated based on rounded up MIC to nearest whole number to match experimental methodology)
SI = (2**test_df['Log2HC10'])/(2**np.ceil(test_df['Log2MIC'])) # Converting Log2 when calculating SI
test_df['SI'] = SI

# Save Total Results dataframe
if reduced_test == 'yes':
    test_df.to_csv(f'{test_path}/Predictions/Reduced Test/Total/Results_{round_num}_Round.csv', encoding='utf-8-sig')
else:
    test_df.to_csv(f'{test_path}/Predictions/All 168000/Total/Results_{round_num}_Round.csv', encoding='utf-8-sig')

TopSI_train_df = train_df.sort_values('SI', ascending=False)[:3]

if bb_train == 'yes':
    
    fig, ax = plt.subplots()
    
    ax.set_box_aspect(1)
    ax.set_xlabel('log2(MIC)')
    ax.set_ylabel('log2(HC10)')
    
    train_bb_df = df[:147]
    train_bb_df = train_bb_df.loc[df['Backbone'] == f'{bb_type}']
    top3SI_train_bb_df = train_bb_df.sort_values('SI', ascending=False)[:3]
    
    ax.scatter(np.log2(train_bb_df['MIC']), np.log2(train_bb_df['HC10']), color = 'blue')
    ax.scatter(np.log2(top3SI_train_bb_df['MIC']), np.log2(top3SI_train_bb_df['HC10']), color = 'orange', edgecolor='black')
    
    ax.plot(log2MIC, log2HC10_SI10, color ='black', alpha=0.3)
    ax.plot(log2MIC, log2HC10_SI20, color ='black', alpha=0.3)
    ax.plot(log2MIC, log2HC10_SI30, color ='black', alpha=0.3)
    
    plt.xticks([1,2,3,4,5,6,7])
    plt.yticks([0,2,4,6,8,10,12,14])

    plt.ylim([0,14.5])
    plt.xlim([0.5,7.5])

    fig.savefig(f'{train_path}/Figures/Top SI Templates/{bb_type}_Templates.svg', format='svg', bbox_inches='tight')
    plt.show()

if design_space_pred == 'yes':
    
    fig, ax = plt.subplots()
    
    ax.set_box_aspect(1)
    ax.set_xlabel('log2(MIC)')
    ax.set_ylabel('log2(HC10)')
    
    ax.scatter(train_df['Log2MIC'], train_df['Log2HC10'], color=[0.7, 0.7, 1.])
    ax.scatter(test_df['Log2MIC'], test_df['Log2HC10'], marker='P', color=[1., 0.7, 0.7])
    
    ax.plot(log2MIC, log2HC10_SI10, color ='black', alpha=0.3)
    ax.plot(log2MIC, log2HC10_SI20, color ='black', alpha=0.3)
    ax.plot(log2MIC, log2HC10_SI30, color ='black', alpha=0.3)
    ax.plot(log2MIC, log2HC10_SI40, color ='black', alpha=0.3)
    ax.plot(log2MIC, log2HC10_SI50, color ='black', alpha=0.3)
    
    ### Plot per round predictions (yellow) vs experimentally validated (green) ###
    
    if round_num == '1st':
        # 1-1
        plt.scatter(4.8366587,7.726872, c='yellow', edgecolors='black', marker='X', label = '1-1')
        plt.scatter(np.log2(128),np.log2(1117.5), c='limegreen', edgecolors='black', marker='X')
        # 1-2
        plt.scatter(5.27304772,8.556483, c='yellow', edgecolors='black', marker='o', label = '1-2')
        plt.scatter(np.log2(256),np.log2(3623), c='limegreen', edgecolors='black', marker='o')
        # 1-3
        plt.scatter(5.36859958,8.613016, c='yellow', edgecolors='black', marker='^', label = '1-3')
        plt.scatter(np.log2(256),np.log2(3457), c='limegreen', edgecolors='black', marker='^')
        
    if round_num == '2nd':
        # 2-1
        plt.scatter(2.626563,4.159363, c='yellow', edgecolors='black', marker='X', label = '2-1')
        plt.scatter(np.log2(32),np.log2(86.5), c='limegreen', edgecolors='black', marker='X')
        # 2-2
        plt.scatter(5.473819,9.828113, c='yellow', edgecolors='black', marker='o', label = '2-2')
        plt.scatter(np.log2(64),np.log2(1481.2), c='limegreen', edgecolors='black', marker='o')
        # 2-3
        plt.scatter(5.364678,9.857255, c='yellow', edgecolors='black', marker='^', label = '2-3')
        plt.scatter(np.log2(64),np.log2(1173.3), c='limegreen', edgecolors='black', marker='^')
        # 2-4
        plt.scatter(4.935137,6.232747, c='yellow', edgecolors='black', marker='d', label = '2-4')
        plt.scatter(np.log2(512),np.log2(4096), c='limegreen', edgecolors='black', marker='d')
        
    if round_num == '3rd':
        # 3-1
        plt.scatter(3.968792,5.66384, c='yellow', edgecolors='black', marker='X', label = '3-1')
        plt.scatter(np.log2(32),np.log2(44.9), c='limegreen', edgecolors='black', marker='X')
        # 3-2
        plt.scatter(3.859946,4.436192, c='yellow', edgecolors='black', marker='o', label = '3-2')
        plt.scatter(np.log2(16),np.log2(13.5), c='limegreen', edgecolors='black', marker='o')
        # 3-3
        plt.scatter(5.713293,10.825716, c='yellow', edgecolors='black', marker='^', label = '3-3')
        plt.scatter(np.log2(256),np.log2(540.3), c='limegreen', edgecolors='black', marker='^')
        # 3-4
        plt.scatter(5.738836,7.464511, c='yellow', edgecolors='black', marker='d', label = '3-4')
        plt.scatter(np.log2(128),np.log2(296.5), c='limegreen', edgecolors='black', marker='d')
        
    if round_num == '4th':
        # 4-1
        plt.scatter(4.967995,9.477851, c='yellow', edgecolors='black', marker='X', label = '4-1')
        plt.scatter(np.log2(128),np.log2(3516.7), c='limegreen', edgecolors='black', marker='X')
        # 4-2
        plt.scatter(4.95616,10.279499, c='yellow', edgecolors='black', marker='o', label = '4-2')
        plt.scatter(np.log2(256),np.log2(4096), c='limegreen', edgecolors='black', marker='o')
        # 4-3
        plt.scatter(4.975095,8.635929, c='yellow', edgecolors='black', marker='^', label = '4-3')
        plt.scatter(np.log2(256),np.log2(3236), c='limegreen', edgecolors='black', marker='^')
        # 4-4
        plt.scatter(5.824251,11.531354, c='yellow', edgecolors='black', marker='d', label = '4-4')
        plt.scatter(np.log2(512),np.log2(4096), c='limegreen', edgecolors='black', marker='d')
        
    if round_num == '5th':
        # 5-1
        plt.scatter(4.967019,8.631038, c='yellow', edgecolors='black', marker='X', label = '5-1')
        plt.scatter(np.log2(32),np.log2(352.1), c='limegreen', edgecolors='black', marker='X')
        # 5-2
        plt.scatter(4.92822,8.841015, c='yellow', edgecolors='black', marker='o', label = '5-2')
        plt.scatter(np.log2(512),np.log2(3902), c='limegreen', edgecolors='black', marker='o')
        # 5-3
        plt.scatter(3.94988,7.564091, c='yellow', edgecolors='black', marker='^', label = '5-3')
        plt.scatter(np.log2(512),np.log2(4096), c='limegreen', edgecolors='black', marker='^')
        # 5-4
        plt.scatter(5.774125,10.234853, c='yellow', edgecolors='black', marker='d', label = '5-4')
        plt.scatter(np.log2(64),np.log2(2178.165942), c='limegreen', edgecolors='black', marker='d')
        
    if round_num == '6th':
        # 6-1
        plt.scatter(6.229602,10.734472, c='yellow', edgecolors='black', marker='X', label = '6-1')
        plt.scatter(np.log2(64),np.log2(2147.83), c='limegreen', edgecolors='black', marker='X')    
        # 6-2
        plt.scatter(6.774514,10.806847, c='yellow', edgecolors='black', marker='o', label = '6-2')
        plt.scatter(np.log2(64),np.log2(1202.26), c='limegreen', edgecolors='black', marker='o')
        # 6-3
        plt.scatter(6.414812,10.949628, c='yellow', edgecolors='black', marker='^', label = '6-3')
        plt.scatter(np.log2(128),np.log2(3556.31), c='limegreen', edgecolors='black', marker='^')
        # 6-4
        plt.scatter(6.462701,10.560432, c='yellow', edgecolors='black', marker='d', label = '6-4')
        plt.scatter(np.log2(64),np.log2(3655.57), c='limegreen', edgecolors='black', marker='d')
    
    plt.xticks([1,2,3,4,5,6,7,8,9])
    plt.yticks([0,2,4,6,8,10,12,14])

    plt.ylim([0,14.5])
    plt.xlim([0.5,10])
    
    plt.legend(loc='lower right', frameon=False, handlelength=0.01, handletextpad=0.6, labelspacing = 0.3, ncol=1)
    
    if reduced_test == 'yes':
        fig.savefig(f'{test_path}/Figures/Design Spaces/Train_Test_{round_num}_Round.svg', format='svg', bbox_inches='tight')

    plt.show()
    
if design_space_sd == 'yes':
    
    plt.rc('font', size=32)          # controls default text sizes
    plt.rc('axes', labelsize=32)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=32)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=32)    # fontsize of the tick labels
    plt.rc('legend', fontsize=32)    # legend fontsize
    plt.rc('figure', titlesize=32)  # fontsize of the figure title
    
    ### HC10 StDev Distribution of Reduced Test Space
    test_df_sd_HC10 = test_df['Log2HC10 s.d.']
    
    fig, ax = plt.subplots(figsize = (11,4))
    
    plt.hist(test_df_sd_HC10, density=True, bins=30, edgecolor='black', color='c')
    plt.xlabel('Log2HC10 Standard Deviation')
    
    ax.set_yticks([])
    ax.set_xticks(np.arange(round(np.min(test_df_sd_HC10),1),round(np.max(test_df_sd_HC10),1),0.2))
    ax.set_xlim(np.min(test_df_sd_HC10),np.max(test_df_sd_HC10))
    ax.xaxis.set_tick_params(width=1.5)
    ax.yaxis.set_tick_params(width=1.5)
    
    for ax, spine in ax.spines.items():
        spine.set_visible(False)
    
    if reduced_test == 'yes':
        fig.savefig(f'{test_path}/Figures/StDev Ranges/HC10_StDev_{round_num}_Round.svg', format='svg', bbox_inches='tight')
    plt.show()
    
    ### MIC StDev Distribution of Reduced Test Space
    test_df_sd_MIC = test_df['Log2MIC s.d.']
    
    fig, ax = plt.subplots(figsize = (11,4))
    
    plt.hist(test_df_sd_MIC, density=True, bins=30, edgecolor='black', color='c')
    plt.xlabel('Log2MIC Standard Deviation')
    
    ax.set_yticks([])
    ax.set_xticks(np.arange(round(np.min(test_df_sd_MIC),1),round(np.max(test_df_sd_MIC),1),0.2))
    ax.set_xlim(np.min(test_df_sd_MIC),np.max(test_df_sd_MIC))
    ax.xaxis.set_tick_params(width=1.5)
    ax.yaxis.set_tick_params(width=1.5)
    
    for ax, spine in ax.spines.items():
        spine.set_visible(False)
    
    if reduced_test == 'yes':
        fig.savefig(f'{test_path}/Figures/StDev Ranges/MIC_StDev_{round_num}_Round.svg', format='svg', bbox_inches='tight')
    plt.show()
        
#####################
##  Round Summary  ##
#####################

# HC10, MIC, Combined Predictions Considered (Table S5)
print(f'\n{len(y_test_HC10)} HC10 predictions')
print(f'{len(y_test_MIC)} MIC predictions')
print((f'{len(test_df)} Total predictions'))

# HC10 and MIC ranges for NSD calculation (Figure 4a)
print(f'\nLog2HC10 StDev Range: {np.round(np.min(test_df["Log2HC10 s.d."]),3)} - {np.round(np.max(test_df["Log2HC10 s.d."]),3)}')
print(f'Log2MIC StDev Range: {np.round(np.min(test_df["Log2MIC s.d."]),3)} - {np.round(np.max(test_df["Log2MIC s.d."]),3)}')
                
                
                
                
                
                