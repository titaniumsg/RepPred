#       Sgourakis Lab
#   Author: Sreeja Kutti Kandy
#   Modified: Sagar Gupta
#   Date: August 30, 2022
#   Email: sreeja@seas.upenn.edu

import pandas as pd
from tqdm import tqdm
from sklearn.utils import shuffle

from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import QuantileTransformer
from sklearn.pipeline import Pipeline
from sklearn.svm import SVR


# import custom libraries
from misc.fasta import FASTA
from misc.constants import *
from misc.common import is_pep_homolog
from create.anchor_class.anchor_class import get_anchor_dict
from create.clustering.clustering_common import read_general_cluster_file

def main():
    
    # Get initial information from (static) database of pHLA structures
    file = FASTA(f"{DATABASE_SCRIPT_PATH}/static/lib_static/PDBS.fasta")
    file.read()
    pdb_dict = file.get_pdb_dict()
    pdb_dict_a0201 = file.get_pdb_dict_allele(selected_allele = "A0201")
    anchor_class = 7; pep_len = 9
    pdb_anchor7_len9_a0201 = get_anchor_dict(pdb_dict_a0201, anchor_class, pep_len)
    a0201_pdbid = list(pdb_anchor7_len9_a0201.keys())

    cluster_info = read_general_cluster_file(f"{DATABASE_SCRIPT_PATH}/static/clustering/greedy_cluster_info_anchor{anchor_class}_len{pep_len}.txt")
    cluster_rep_pdbid = []
    for num_pdbid in cluster_info.keys():
        cluster_rep_pdbid.append(num_pdbid.split("_")[1])
    
    # Read energy and distance score data from threading benchmark
    data = pd.read_csv('score_parameters_withdxtal.csv')
    data = data.set_index(data.columns[0], drop=True)
    data = data.loc[:, data.var() > 0.00]
    data.rename(columns={"Unnamed: 132": "dscore_xtal"}, inplace=True)

    true_column_name = 'dscore_xtal'

    a0201_rep_index = []                               
    for ind in data.index:
        target = ind.split('_on_')[0]
        template = ind.split('_on_')[1]
        if (target in a0201_pdbid) and (template in cluster_rep_pdbid) and (target not in cluster_rep_pdbid): # Filter data to keep only A0201 x representatives
            target_pep = pdb_dict[target].split("_")[0]
            template_pep = pdb_dict[template].split("_")[0]
            if not is_pep_homolog(target_pep, template_pep)[0]:  # remove target x representatives pairs that are homologs
                a0201_rep_index.append(ind)

    data_a0201_rep = data.reindex(a0201_rep_index)

    a0201_targets = list(set([pair[0] for pair in data_a0201_rep.index.str.split('_on_', expand=True)]))

    print(f"There are {len(a0201_pdbid)} A0201 structures in the database and our benchmark has {len(a0201_targets)} A0201 structures")
    print(f"There are {len(cluster_rep_pdbid)} representative structures")

    with open("a0201_targets.txt", "w") as txtfile:
        for pdbid in a0201_targets:
            txtfile.write(pdbid+"\n")

    def leave_one_out_split(data, target_id):

        train_index = []
        test_index = []

        for ind in data.index:
            target = ind.split('_on_')[0]
            template = ind.split('_on_')[1]
            if target != target_id and template != target_id and target != template:
                train_index.append(ind)
            
            elif target == target_id and template != target_id:
                test_index.append(ind)

        train_data = data.reindex(train_index)
        test_data = data.reindex(test_index)

        return train_data, test_data

    y_pred_store = pd.DataFrame()
    
    for target_id in tqdm(a0201_targets):

        # Regression model information
        transformer = QuantileTransformer(output_distribution='uniform')
        svr = SVR(kernel='rbf', cache_size=2000, C=10, epsilon=0.1)
        qsvr = Pipeline(steps=[('t', transformer), ('m', svr)])
        
        # Get training data
        train_data, test_data = leave_one_out_split(data_a0201_rep, target_id)
        train_data = shuffle(train_data)
        X_train = train_data.drop(['totalScore', 'peptideScore', true_column_name], axis=1)
        y_train = train_data[true_column_name]

        # Get test data
        X_test = test_data.drop(['totalScore','peptideScore', true_column_name],axis=1)
        scaler = StandardScaler().fit(X_train)
        X_train_scaled = scaler.transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # Predict D-score values
        qsvr.fit(X_train_scaled, y_train)
        y_pred = qsvr.predict(X_test_scaled)
        test_data['dscore_pred'] = y_pred
        
        # Store prdicted D-score values
        y_pred_store = pd.concat([y_pred_store, test_data[[true_column_name, 'dscore_pred']][:]], ignore_index=False)

    y_pred_store.to_csv("y_pred_nohomolog.csv", index=True)

if __name__ == "__main__":

    main()


