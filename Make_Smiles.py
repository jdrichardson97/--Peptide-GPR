# -*- coding: utf-8 -*-
"""
Custom script to turn a/B-peptide sequences
into their SMILES string equivalent
(see Figure S2)

Author: Joshua D Richardson
"""

import numpy as np
import pandas as pd
import os

#############################################################################################################################
# BEGIN USER INPUT PARAMETERS

# Define these according to your local directory #
user = 'jdrichardso3'
home_dir = f'C:/Users/{user}/OneDrive - UW-Madison/Documents/Research/Papers/Paper Drafts/2_aB Peptide Predictions/Data/'

# Input parameters
bb = 'aaB' # aaB (Rounds 1-4) or aaBaaaB (Rounds 5-6)

# END USER INPUT PARAMETERS
##############################################################################################################################

# Training and Test Paths based on home_dir above
train_path = f'{home_dir}/Training'
test_path = f'{home_dir}/Test Design Space'

# 20 possible alpha and 21 possible beta residues to vary on templates (see Figure 1c-d)
alpha_res = ['Gly','Ile','Leu','Ala','Lys','Arg','Asp','Glu','Ser','Thr','Phe','Val','Abu','Met','Trp','Tyr','Asn','Gln','Nva','Nle']
beta_res = ['βGly','βIle','βLeu','βAla','βLys','βArg','βAsp','βGlu','βSer','βThr','βPhe','βVal','βAbu','βMet','βTrp','βTyr','βAsn','βGln','βNva','βNle', 'ACPC']

# All peptide sidechain fragments (see Figure S2b)
side = {# Present in Expanded Dataset
        'Gly': '',
        'Ile': '[C@H](CC)C',
        'Leu': 'CC(C)C',
        'Ala': 'C',
        'Lys': 'CCCC[NH3+]',
        'Arg': 'CCCNC(=[NH2+])N',
        'Asp': 'CC(=O)[O-]', #CC(=O)O
        'Glu': 'CCC(=O)[O-]', #CCC(=O)O
        'Ser': 'CO',
        'Thr': '[C@H](C)O',
        'Phe': 'Cc9ccccc9',
        'Val': 'C(C)C',
        'Abu': 'CC',
        'ACPC': 'CCC', # Need to add number at end (CCC#) depending on place in backbone (see below)
        # Not Present (Natural) 
        'Met':'CCSC',
        'Trp':'Cc8c[nH]c9c8cccc9',
        #'Pro':'',
        #'Cys':'',
        'Tyr':'Cc8ccc(cc8)O',
        'Asn':'CC(=O)N',
        'Gln':'CCC(=O)N',
        #'His':'',
        # Not Present (Synthetic)
        'Nva':'CCC',
        'Nle':'CCCC',
        }

seq = []
smiles = []

# Templating new sequences for Rounds 1-4 on aaB backbone (Figure 1c)
if bb == 'aaB':
    
    Res123_pos = alpha_res
    Res4_pos = beta_res
    
    if os.path.exists(f'{test_path}/SMILES_aaB.csv') == False:
        # Keeping track of iteration progress
        tot_it = (len(Res123_pos)**3)*len(Res4_pos)
        it = 0
        
        for a in Res123_pos:
            for b in Res123_pos:
                for c in Res123_pos:
                    for d in Res4_pos:
                        pep = 'RES1-Leu-ACPC-Lys-RES2-β3hAla-Lys-Lys-ACPC-RES3-Lys-RES4-Phe-CONH2'
                        pep = pep.replace('RES1',a)
                        pep = pep.replace('RES2',b)
                        pep = pep.replace('RES3',c)
                        pep = pep.replace('RES4',d)

                        seq = np.append(seq, pep)
                        
                        pep_list = pep.split('-')
                        pep_list = pep_list[:-1]
                        smi_start = '[NH3+]'
                        beta = [] # Fill with 0 if not beta and 1 if beta for standard amino acid residues
                        aa = [] # aa identities for side chains
                        
                        # Append to beta and aa lists
                        for i in pep_list:
                            #j = '1'
                            if i[0:3] == 'β':
                                beta = np.append(beta,1)
                                aa = np.append(aa, i[3:])
                            else:
                                beta = np.append(beta,0)
                                aa = np.append(aa, i)
                        
                        # Create Backbone Portion of Smiles String
                        y = 1
                        for x in range(len(aa)):
                            if aa[x] == 'ACPC':
                                smi_start += f'[C@H]{y}[C@H](C(=O)N'
                                y += 1
                            elif aa[x] == 'Gly' and beta[x] == 0:
                                smi_start += 'CC(=O)N'
                            elif aa[x] == 'Gly' and beta[x] == 1:
                                smi_start += 'CCC(=O)N'
                            elif aa[x] != 'Gly' and beta[x] == 0:
                                smi_start += '[C@H](C(=O)N'
                            elif aa[x] != 'Gly' and beta[x] == 1:
                                smi_start += '[C@H](CC(=O)N'
                            elif aa[x] == 'ACPC':
                                smi_start += f'[C@H]{y}[C@H](C(=O)N'
                                
                        smi_bb = smi_start + ')'
                        
                        # Add Sidechains
                        y-=1
                        smi_pep = smi_bb
                        aa_rev = np.flip(aa)
                        for z in range(len(aa_rev)):
                            if z != len(aa_rev):
                                if aa_rev[z] == 'Gly':
                                    smi_pep = smi_pep + side[aa_rev[z]]
                                elif aa_rev[z] == 'ACPC':
                                    smi_pep = smi_pep + side[aa_rev[z]] + str(y) + ')'
                                    y-=1
                                else:
                                    smi_pep = smi_pep + side[aa_rev[z]] + ')'
                            else:
                                if aa_rev[z] != 'ACPC':
                                    smi_pep = smi_pep + side[aa_rev[z]]
                                else:
                                    smi_pep = smi_pep + side[aa_rev[z]] + str(y)
                                    y-=1
                        
                        # Remove extraneous parenthesis at end (mainly caused by Gly being first residue)
                        while smi_pep[-1] == ')':
                            smi_pep = smi_pep.removesuffix(')')
                        
                        # Final smi string
                        smiles = np.append(smiles, smi_pep)
                        it+=1
                        if it%1000 == 0:
                            print(f'Done with iteration {it} of {tot_it}')
    
        #Create Smiles Dataframe for export
        smi_df = pd.DataFrame(index=np.arange(tot_it), columns=['seq', 'smi'])
        smi_df['seq'] = seq
        smi_df['smi'] = smiles
        
        smi_df.to_csv(f'{test_path}/SMILES_aaB.csv')
    
    else:
        smi_df = pd.read_csv(f'{test_path}/SMILES_aaB.csv', index_col = 0, header = 0)
        seq = smi_df['seq']
        smiles = smi_df['smi']
    
# Templating new sequences for Rounds 5-6 on aaBaaaB backbone (Figure 1d)
elif bb == 'aaBaaaB': 
    
    Res134_pos = alpha_res
    Res2_pos = beta_res
    
    if os.path.exists(f'{test_path}/SMILES_aaBaaaB.csv') == False:
        # Keeping track of iteration progress
        tot_it = (len(Res134_pos)**3)*len(Res2_pos)
        it = 0
        
        for a in Res134_pos:
            for b in Res2_pos:
                for c in Res134_pos:
                    for d in Res134_pos:
                        pep = 'RES1-ACPC-Phe-Lys-Ile-RES2-Lys-Lys-ACPC-RES3-Lys-RES4-βPhe-CONH2'
                        pep = pep.replace('RES1',a)
                        pep = pep.replace('RES2',b)
                        pep = pep.replace('RES3',c)
                        pep = pep.replace('RES4',d)

                        seq = np.append(seq, pep)
                        
                        pep_list = pep.split('-')
                        pep_list = pep_list[:-1]
                        smi_start = '[NH3+]'
                        beta = [] # Fill with 0 if not beta and 1 if beta for standard amino acid residues
                        aa = [] # aa identities for side chains
                        
                        # Append to beta and aa lists
                        for i in pep_list:
                            #j = '1'
                            if i[0:3] == 'β':
                                beta = np.append(beta,1)
                                aa = np.append(aa, i[3:])
                            else:
                                beta = np.append(beta,0)
                                aa = np.append(aa, i)
                        
                        # Create Backbone Portion of Smiles String
                        y = 1
                        for x in range(len(aa)):
                            if aa[x] == 'ACPC':
                                smi_start += f'[C@H]{y}[C@H](C(=O)N'
                                y += 1
                            elif aa[x] == 'Gly' and beta[x] == 0:
                                smi_start += 'CC(=O)N'
                            elif aa[x] == 'Gly' and beta[x] == 1:
                                smi_start += 'CCC(=O)N'
                            elif aa[x] != 'Gly' and beta[x] == 0:
                                smi_start += '[C@H](C(=O)N'
                            elif aa[x] != 'Gly' and beta[x] == 1:
                                smi_start += '[C@H](CC(=O)N'
                            elif aa[x] == 'ACPC':
                                smi_start += f'[C@H]{y}[C@H](C(=O)N'
                                
                        smi_bb = smi_start + ')'
                        
                        # Add Sidechains
                        y-=1
                        smi_pep = smi_bb
                        aa_rev = np.flip(aa)
                        for z in range(len(aa_rev)):
                            if z != len(aa_rev):
                                if aa_rev[z] == 'Gly':
                                    smi_pep = smi_pep + side[aa_rev[z]]
                                elif aa_rev[z] == 'ACPC':
                                    smi_pep = smi_pep + side[aa_rev[z]] + str(y) + ')'
                                    y-=1
                                else:
                                    smi_pep = smi_pep + side[aa_rev[z]] + ')'
                            else:
                                if aa_rev[z] != 'ACPC':
                                    smi_pep = smi_pep + side[aa_rev[z]]
                                else:
                                    smi_pep = smi_pep + side[aa_rev[z]] + str(y)
                                    y-=1
                        
                        # Remove extraneous parenthesis at end (mainly caused by Gly being first residue)
                        while smi_pep[-1] == ')':
                            smi_pep = smi_pep.removesuffix(')')
                        
                        # Final smi string
                        smiles = np.append(smiles, smi_pep)
                        it+=1
                        if it%1000 == 0:
                            print(f'Done with iteration {it} of {tot_it}')
    
        #Create Smiles Dataframe for export
        smi_df = pd.DataFrame(index=np.arange(tot_it), columns=['seq', 'smi'])
        smi_df['seq'] = seq
        smi_df['smi'] = smiles
        
        smi_df.to_csv(f'{test_path}/SMILES_aaBaaaB.csv')
    
    else:
        smi_df = pd.read_csv(f'{test_path}/SMILES_aaBaaaB.csv', index_col = 0, header = 0)
        seq = smi_df['seq']
        smiles = smi_df['smi']

   