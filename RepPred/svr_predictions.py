#!/usr/bin/env python

import subprocess
import pandas as pd
import pickle as pkl
import shutil
import os

os.system("rm scores/model_scores.csv")

combine=subprocess.run(f"bash combine_scores.sh",  shell=True)

def get_target_template_pairs(pep,datax):              # Get data for train and test 
    keepindex = []
    for ind in datax.index:
        indsp = ind.split('_on_')
        if (indsp[0]==pep): 
            keepindex += [ind]
    return(datax.reindex(keepindex))

data_p=pd.read_csv('scores/model_scores.csv',sep=',')
data_p=data_p.set_index(data_p.columns[0], drop=True)
peptide_list=[]
for ind in data_p.index:
    indsp = ind.split('_on_')
    peptide_list.append(indsp[0])
peptide_list=list(set(peptide_list))
print(peptide_list)

rmodel = pkl.load(open('reg_model.pkl', 'rb'))
rscaler = pkl.load(open('reg_scaler.pkl', 'rb'))

predictions=pd.DataFrame()
for n,peptide in enumerate(peptide_list):
    pep_models=get_target_template_pairs(peptide,data_p)
    X_test=pep_models.drop(['totalScore','peptideScore'],axis=1)
    X_test_scaled = rscaler.transform(X_test)
    pep_models['ypred']=rmodel.predict(X_test_scaled)
    pep_models_sorted=pep_models.sort_values(by=['ypred'])
    pep_models_score_sorted=pep_models.sort_values(by=['totalScore'])

    predictions.loc[n,'target']=peptide
    predictions.loc[n,'sequence']=open(f'./input_sequence/{peptide}_seq.txt','r').readlines()[1].strip()
    predictions.loc[n,'best_dscore']=pep_models_sorted['ypred'][0]
    predictions.loc[n,'best_model']=pep_models_sorted.index[0]
    shutil.copyfile(f'run-dir/relaxed_structures/best_score/{pep_models_sorted.index[0]}.pdb', f'model_predictions/{peptide}_model.pdb',follow_symlinks=True)
predictions.to_csv(f'model_predictions/predictions_{peptide}.csv', index=True, mode='a')

