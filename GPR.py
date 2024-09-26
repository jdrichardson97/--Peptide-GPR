# -*- coding: utf-8 -*-
"""
Script to run GPR model training and test design 
space HC10/MIC predictions

Author: Joshua D Richardson
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold, GridSearchCV, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Lasso, LassoCV
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import DotProduct, RBF, RationalQuadratic, ExpSineSquared

# Set graphing properties
plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize
plt.rc('figure', titlesize=18)  # fontsize of the figure title

#############################################################################################################################
# BEGIN USER INPUT PARAMETERS

# Define these according to your local directory #
user = 'jdrichardso3'
home_dir = f'C:/Users/{user}/OneDrive - UW-Madison/Documents/Research/Papers/Paper Drafts/2_aB Peptide Predictions/Data/'

# Input parameters
label = 'HC10' # HC10 MIC
round_num = '6th' # 1st 2nd 3rd 4th 5th 6th
update_desc = 'yes' # Whether to update descriptors (1-3 = no, 4-6 = yes)
eval_on_test = 'yes' # Whether to evaluate trained model on 168,000 sequence test design space
desc_red = 'yes' # Reduce test descriptors to be in range of training descriptors (default is yes, see Fig S7-S8)
cv_num = 10 # Number of folds for cross-validation (CV) --> 10 folds were used in this study (see Figure S5)

# GPR Input Parameters
select_desc = 'auto' # Default is 'auto' to allow Gridsearch to select hyperparams (Table S4)

if select_desc == 'manual': # Pick hyperparams manually (e.g., evaluate model robustness --> see Figure S12)
    gpr_kappa = 0.01
    gpr_kernel = ExpSineSquared(length_scale=0.01, periodicity=100)
    
# END USER INPUT PARAMETERS
##############################################################################################################################

# Training and Test Paths based on home_dir above
train_path = f'{home_dir}/Training'
test_path = f'{home_dir}/Test Design Space'

print(f'\nWorking on {label}, round = {round_num}, update_desc = {update_desc}\n')

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

# Import training peptides (Figure 2a)
df = pd.read_excel(f'{train_path}/Sequence_Activity.xlsx', header=0)
pep_num = df['Index'].dropna()
seq = df['Sequence'][:len(pep_num)]

# All 200 RDKit training descriptors and labels
data_raw = pd.read_csv(
    f'{train_path}/RDKit.csv', header=0)
data = data_raw
data = data.head(row_keep) # Number of sequences available based on prediction round

# Remove Large RDKit Ipc descriptor causing problems with Pearson/Regression
if 'Ipc' in data.columns:
    data = data.drop(labels='Ipc', axis=1)
    
# Remove constant descriptors (Figure S3)
data = data.loc[:, (data != data.iloc[0]).any()]

# Import test sequence space peptides
if round_num=='1st' or round_num=='2nd' or round_num=='3rd' or round_num=='4th':
    new_df = pd.read_csv(
        f'{test_path}/SMILES_aaB.csv', header=0, index_col=0)  # aaB backbone for Rounds 1-4
else:
    new_df = pd.read_csv(
        f'{test_path}/SMILES_aaBaaaB.csv', header=0, index_col=0)  # aaBaaaB backbone for Rounds 5-6
new_seq = new_df['Sequence']
            
# Proportionate Allocation of Labels (Figure S5d)
def df_propalloc(data, n_splits):
    set_size = int(len(data.index)/n_splits)
    gap = n_splits
    
    num2sort = len(data)
    tot = set_size * gap
    remain = num2sort - tot
    remain_cols = data.tail(remain)
    
    new_df = pd.DataFrame(columns = data.columns)
    for i in range(n_splits):
        for j in range(set_size):
            new_df = new_df.append(data.loc[[i]])
            i += gap
    
    def Insert_row_(row_number, df, row_value): # Function to insert remaining rows into dataframe

        df1 = df[0:row_number]
        df2 = df[row_number:]
        df1.loc[row_number] = row_value
      
        df_result = pd.concat([df1, df2])
  
        return df_result
    
    df_result = new_df
    df_result = df_result.reset_index(drop=True)
    if remain != 0:
        insert_num = set_size
        add = 1
        mult = 2
        for i in range(len(remain_cols)):
            df_result = Insert_row_(insert_num, df_result, remain_cols.iloc[i])
            df_result = df_result.reset_index(drop=True)
            insert_num = set_size*mult + add
            mult +=1
            add +=1
        new_df = df_result
    
    new_df = new_df.reset_index(drop=True) #Reset indices if needed
    return new_df

sorted_df = data.sort_values(label)
sorted_df = sorted_df.reset_index(drop=True)
n_splits = cv_num # Define for CV
data = df_propalloc(sorted_df, n_splits)
    
# Defining label properties
y_df = data.loc[:, data.columns == label]
y = np.ravel(y_df)
    
# Log 2 Scale label (Figure S6)
y_df = y_df.astype(float)
y_df = np.log2(y_df)
y = y.astype(float)
y = np.log2(y)
scale = 'Log2'

# Remove remaining non-descriptor columns from RDKit descriptor data and det descriptor matrix (X)
if 'Index' in data.columns:
    index_df = data.loc[:, data.columns == 'Index']
    index = np.ravel(index_df)
    data = data.drop(labels = 'Index', axis=1)
if 'MIC' in data.columns:
    mic_df = data.loc[:, data.columns == 'MIC']
    mic = np.ravel(mic_df)
    data = data.drop(labels = 'MIC', axis=1)  
if 'HC10' in data.columns:
    hc10_df = data.loc[:, data.columns == 'HC10']
    hc10 = np.ravel(hc10_df)
    data = data.drop(labels = 'HC10', axis=1)

X_df = data
X = X_df.to_numpy()

# Whether to update descriptors for Prediction Round (Figure S4)
if update_desc == 'yes':
    # Assign standardized matrices for Lasso CV
    scaler = StandardScaler()
    X_LASSO = scaler.fit_transform(X)
    X_LASSO_df = pd.DataFrame(scaler.fit_transform(X_df), columns=X_df.columns)
    
    # Set cross validation procedure and do LassoCV to remove unimportant descriptors
    alpha_range = np.logspace(-3,0,100)
    
    # Edge case for 6th round MIC
    if label == 'MIC':
        if round_num == '6th':
            alpha_range = alpha_range[3:]
    
    cv = KFold(n_splits=cv_num, shuffle=False)
    lassocv = LassoCV(cv=cv, alphas = alpha_range, max_iter = 1000000)
    lassocv.fit(X_LASSO, y)
    
    # Splits for Lasso CV
    for i, (train_index, test_index) in enumerate(cv.split(X_LASSO)):
        print(f"Fold {i}:")
        #print(f"  Train: index={train_index}")
        print(f"  Test:  index={test_index}")
        print(f"  Test:  values={y[test_index]}")
    
    # Lasso using alpha determined by Lasso CV
    
    lasso = Lasso(max_iter=1000000)
    lasso.set_params(alpha=lassocv.alpha_)
    print("Alpha=", lassocv.alpha_)
    lasso.fit(X_LASSO, y)
    print("Top 5 Kept coefficients:")
    Coeffs = pd.Series(lasso.coef_, index=X_df.columns) #
    
    opt_alpha = lassocv.alpha_
    
    #Set close to zero coefficients as zero
    for i in range(len(Coeffs)):
        if abs(Coeffs[i]) < 10**-5: #10^-5
            Coeffs[i] = 0
            
    # Plotting MSE vs Alpha
    fig, ax = plt.subplots()
    ax.set_box_aspect(1)
    #ax.set_title("LASSO Cross Validation")
    ax.set_xlabel("α")
    ax.set_ylabel("RMSE")
    rmse_path = np.sqrt(lassocv.mse_path_)
    ax.semilogx(lassocv.alphas_, rmse_path, ":", linewidth=2)
    ax.plot(
        lassocv.alphas_ ,
        rmse_path.mean(axis=-1),
        "k",
        label="Average",
        linewidth=2.5)
    ax.axvline(lassocv.alpha_, linestyle="--", color="k", label=f"α = {np.round(opt_alpha,4)}", linewidth=2.5)
    ax.xaxis.set_tick_params(width=1.5)
    ax.yaxis.set_tick_params(width=1.5)
    ax.axis("tight")
    plt.xticks([10**-3,10**-2,10**-1,1])
    fig.savefig(f'{train_path}/Figures/LASSO CV/{label}_{round_num}_Round.svg', format='svg', bbox_inches='tight')
    
    plt.show()
    
    Kept_Coeffs_pre = abs(Coeffs.iloc[Coeffs.to_numpy().nonzero()])
    Kept_Coeffs = Kept_Coeffs_pre.sort_values(ascending=False) # Sort from most to least important
    
    print(Kept_Coeffs[:5]) # Top 5 Coefficients kept
    print('Kept: ', np.shape(Kept_Coeffs)[0])
    print('Removed: ', (np.shape(Coeffs)[0] - np.shape(Kept_Coeffs)[0]))
    Kept_Coeffs_HC10 = Kept_Coeffs
    
    Kept_Coeffs.to_csv(f'{train_path}/LASSO CV Coeffs/Kept_{label}_Coeffs_{round_num}_Round.csv', encoding='utf-8-sig')

    # Removing unimportant descriptors before modelling (after Lasso CV)
    remove_coeff = list(set(X_df.columns).difference(Kept_Coeffs.index.values))
    X_df = X_df[list(Kept_Coeffs.index)]
    X = X_df.to_numpy()

elif update_desc == 'no':
    Kept_Coeffs = pd.read_csv(f'{train_path}/LASSO CV Coeffs/Kept_{label}_Coeffs_1st_Round.csv', header=0)
    X_df = X_df[list(Kept_Coeffs['Unnamed: 0'])]
    X = X_df.to_numpy()

X_df_train = X_df

# Import new test sequence rdkit descriptors
if round_num == '1st' or round_num == '2nd' or round_num == '3rd' or round_num == '4th':
    data_raw_test = pd.read_csv(f'{test_path}/RDKit_aaB.csv', header=0)
else:  
    data_raw_test = pd.read_csv(f'{test_path}/RDKit_aaBaaaB.csv', header=0)

# Match test descriptors with training descriptors
X_df_test_pre = data_raw_test.loc[:, data_raw_test.columns != 'Index']
X_df_test_pre = X_df_test_pre[X_df.columns]
X_df_test = X_df_test_pre

# Descriptor reduction within range of training (Figure S7-S8)
if desc_red == 'yes':
    
    if update_desc == 'no':
                
        if os.path.exists(f'{test_path}/Reduced Test Spaces/{label}/inrange_for_{round_num}_Round_orig.csv') == False:
            def minMax(x):
                return pd.Series(index=['min','max'],data=[x.min(),x.max()])
        
            train_minmax = X_df_train.apply(minMax)
            test_bool = X_df_test.copy()
            
            idx=1
            for i in X_df_test.columns:
                for j in range(len(X_df_test.index)):
                    a =  (train_minmax[i].iloc[0] <= test_bool[str(i)].iloc[j] <= train_minmax[i].iloc[1])
                    test_bool[i].iloc[j] = a
                print(f'Done with test matrix column {idx} of {len(X_df_test.columns)}')
                idx+=1
                
            result = test_bool.all(axis='columns')
            inrange_idx = list(result[result].index.values)
            inrange_idx_df = pd.DataFrame(inrange_idx)
            inrange_idx_df.to_csv(f'{test_path}/Reduced Test Spaces/{label}/inrange_for_{round_num}_Round_orig.csv', encoding='utf-8-sig')
            X_df_test = X_df_test.loc[list(inrange_idx)]
            
        else:
            inrange_idx_df = pd.read_csv(f'{test_path}/Reduced Test Spaces/{label}/inrange_for_{round_num}_Round_orig.csv', header=0)
            inrange_idx = inrange_idx_df['0'].values.flatten()
            X_df_test = X_df_test.loc[list(inrange_idx)]
                
    elif update_desc == 'yes':
        
        if os.path.exists(f'{test_path}/Reduced Test Spaces/{label}/inrange_for_{round_num}_Round_update.csv') == False:
            def minMax(x):
                return pd.Series(index=['min','max'],data=[x.min(),x.max()])
        
            train_minmax = X_df_train.apply(minMax)
            test_bool = X_df_test.copy()
            
            idx=1
            for i in X_df_test.columns:
                for j in range(len(X_df_test.index)):
                    a =  (train_minmax[i].iloc[0] <= test_bool[str(i)].iloc[j] <= train_minmax[i].iloc[1])
                    test_bool[i].iloc[j] = a
                print(f'Done with test matrix column {idx} of {len(X_df_test.columns)}')
                idx+=1
                
            result = test_bool.all(axis='columns')
            inrange_idx = list(result[result].index.values)
            inrange_idx_df = pd.DataFrame(inrange_idx)
            inrange_idx_df.to_csv(f'{test_path}/Reduced Test Spaces/{label}/inrange_for_{round_num}_Round_update.csv', encoding='utf-8-sig')
            X_df_test = X_df_test.loc[list(inrange_idx)]
            
        else:
            inrange_idx_df = pd.read_csv(f'{test_path}/Reduced Test Spaces/{label}/inrange_for_{round_num}_Round_update.csv', header=0)
            inrange_idx = inrange_idx_df['0'].values.flatten()
            X_df_test = X_df_test.loc[list(inrange_idx)]

# Standardize Training Matrix
scaler = StandardScaler()

X_df_train = pd.DataFrame(scaler.fit_transform(X_df_train), columns=X_df_train.columns)
y_df_train = y_df
y_df_train = y_df_train.reset_index(drop=True)
y_train = y_df_train.to_numpy().flatten()

# Apply scaling to test set (if applicable)
if eval_on_test == 'yes':
    X_df_test = pd.DataFrame(scaler.transform(X_df_test), columns=X_df_test.columns)
    if desc_red == 'yes':
        X_df_test = X_df_test.set_index([pd.Index(inrange_idx)])
    elif desc_red == 'no':
        X_df_test = X_df_test.set_index([pd.Index([i for i in range(0,168000)])])

# GPR Gridsearch
if select_desc == 'auto':
    ##Gaussian Process Regressor
    gaus_model = GaussianProcessRegressor(optimizer='fmin_l_bfgs_b')
    
    gaus_vals = [i for i in np.logspace(-2,2,5)]
    gaus_alpha = [1e-20, 1e-15, 1e-10, 1e-5, 0.01, 0.1, 1.0, 5.0, 10.0]
    
    gaus_param = [{
        "alpha":  gaus_alpha,
        "kernel": [RBF(length_scale=k) for k in gaus_vals]},
    {
        "alpha":  gaus_alpha,
        "kernel": [RationalQuadratic(length_scale=k, alpha=l) for k in gaus_vals for l in gaus_vals]},
    {
        "alpha":  gaus_alpha,
        "kernel": [ExpSineSquared(length_scale=k, periodicity=l) for k in gaus_vals for l in gaus_vals]},
    {
        "alpha":  gaus_alpha,
        "kernel": [DotProduct(sigma_0=k) for k in gaus_vals]}
        ]
    
    grid_gaus = GridSearchCV(estimator=gaus_model, param_grid=gaus_param ,n_jobs=-1, cv=cv_num)
    grid_gaus.fit(X_df_train, y_df_train)

    gaus = GaussianProcessRegressor(
                        alpha = grid_gaus.best_params_['alpha'],
                        kernel = grid_gaus.best_params_['kernel'],
                        )
    
elif select_desc == 'manual':
    gaus = GaussianProcessRegressor(
                        alpha = gpr_kappa,
                        kernel = gpr_kernel,
                        )
    
gaus.fit(X_df_train, y_df_train)
y_predicted_gaus = cross_val_predict(gaus, X_df_train, y_df_train, cv=cv_num)
y_predicted_gaus = y_predicted_gaus.flatten()

# Parity plot visualization on training

# Recover original ordering of sequences (if proportionately allocated at the beginning)
idx_list = [int(i) for i in index]
y_train_reindex_df = pd.DataFrame(y_train, index=idx_list).sort_index()
y_test_reindex_df = pd.DataFrame(y_predicted_gaus, index=idx_list).sort_index()
y_train_ri = y_train_reindex_df.values
y_test_ri = y_test_reindex_df.values

sf = 0.3
fig, ax = plt.subplots()
ax.scatter(y_train_ri[:147], y_test_ri[:147], label='R\u00b2 = {:.3f}'.format(r2_score(y_train, y_predicted_gaus)), s=40, facecolors='none', edgecolors='blue')
ax.scatter(y_train_ri[147:], y_test_ri[147:], s=50, c='red', marker='^')

ax.set_xlabel(f'Actual Log2{label}')
ax.set_ylabel(f'Predicted Log2{label}')

#Set Parity Line
if label == 'HC10':
    ax.plot([0, 12.5], [0, 12.5], "k--", lw=4)
if label == 'MIC':
    ax.plot([0, 10], [0, 10], "k--", lw=4)

ax.set_box_aspect(1)

ax.legend(frameon=False, loc='lower right', handlelength=0.3, handletextpad=0.6, labelspacing = 0.3)

if label == 'HC10':
    plt.xlim([0, 12.5])
    plt.ylim([0, 12.5])
    plt.xticks(range(0,13,2))
    plt.yticks(range(0,13,2))
    ax.xaxis.set_tick_params(width=1.5)
    ax.yaxis.set_tick_params(width=1.5)
    
elif label == 'MIC':
    plt.xlim([1, 10])
    plt.ylim([1, 10])
    plt.xticks(range(2,10,1))
    plt.yticks(range(2,10,1))
    ax.xaxis.set_tick_params(width=1.5)
    ax.yaxis.set_tick_params(width=1.5)

fig.savefig(f'{train_path}/Figures/Parity Plots/{label}_{round_num}_Round.svg', format='svg', bbox_inches='tight')
plt.show()
    
# Output y_train and y_test to pandas df
train_test_df = pd.DataFrame([])
train_test_df.index += 1
train_test_df['Train'] = y_train_ri.flatten()
train_test_df['Predicted'] = y_test_ri.flatten()
train_test_df.to_csv(f'{train_path}/CV Predict/{label}_{round_num}_CV_predict.csv', encoding='utf-8-sig')

if eval_on_test == 'yes':
    # Predict new seq HC10 or MIC
    if select_desc == 'auto':
        gaussian = GaussianProcessRegressor(
                            alpha = grid_gaus.best_params_['alpha'],
                            kernel = grid_gaus.best_params_['kernel'],
                          )
    elif select_desc == 'manual':
        gaussian = GaussianProcessRegressor(
                            alpha = gpr_kappa,
                            kernel = gpr_kernel,
                          )
    gaussian.fit(X_df_train,y_df_train)
    
    y_pred, y_std = gaussian.predict(X_df_test, 
                                   return_std=True
                                   )
    
    # Final Results Df
    results_df = pd.DataFrame(index = X_df_test.index)
    results_df['Sequence'] = new_seq.iloc[list(X_df_test.index)]
    results_df[label] = y_pred
    results_df['Stdev'] = y_std
    
    if desc_red == 'yes':
        
        results_df.to_csv(f'{test_path}/Predictions/Reduced Test/{label}/{label}_Results_{round_num}_Round.csv', encoding='utf-8-sig')
        
        if round_num == '4th' and update_desc == 'no': # Hypothetical Round 4 with same kept descriptors as Round 1-3 (used in Fig S14)
        
            results_df.to_csv(f'{test_path}/Predictions/Reduced Test/{label}/{label}_Results_{round_num}_Round_noupdate.csv', encoding='utf-8-sig')
            
    elif desc_red == 'no':
            
        results_df.to_csv(f'{test_path}/Predictions/All 168000/{label}/{label}_Results_{round_num}_Round.csv', encoding='utf-8-sig')
          
#####################
##  Round Summary  ##
#####################

# Descriptors used (Table S1 or S2)
print('\n')
print('Descriptors:\n')
print(Kept_Coeffs)

# Trained GPR model parameters (Table S4)
print('\n')
print('GPR Parameters: ')
if select_desc == 'auto':
    print(f"Κ = {grid_gaus.best_params_['alpha']}")
    print(f"kernel = {grid_gaus.best_params_['kernel']}")
elif select_desc == 'manual':
    print(f"Κ = {gpr_kappa}")
    print(f"kernel = {gpr_kernel}")
    
# Test Sequence Design space (Table S5)
if eval_on_test == 'yes':
    print('\n')
    print(f'Test Sequences Considered: {len(results_df)}')

