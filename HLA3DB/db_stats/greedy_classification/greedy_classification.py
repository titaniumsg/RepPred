#       Sgourakis Lab
#   Author: Sagar Gupta
#   Date: May 14, 2022
#   Email: sagarg@sas.upenn.edu

# import required libraries
from collections import defaultdict

# import custom libraries
from misc.fasta import FASTA
from misc.constants import *
from misc.common import *

def get_binary_matrix(pdb_dihedrals, anchor_class):

    binary_matrix = []

    for pdb1 in pdb_dihedrals.keys():
        row = []
        for pdb2 in pdb_dihedrals.keys():
            num_angles = anchor_class*2-6

            if distance_score_threshold2(pdb1, pdb2, pdb_dihedrals, SUM_DIST_THRESHOLD/8*num_angles)[0]:
                row.append(1)
            else:
                row.append(0)

        binary_matrix.append(row)

    return binary_matrix

def get_greedy_class(binary_matrix, pdbids):

    sum_values = {}

    for index, row in enumerate(binary_matrix):
        sum_values[pdbids[index]] = sum(row)

    rep = max(sum_values, key=sum_values.get)
    rep_index = pdbids.index(rep)

    neighbors, neighbors_index = [], []

    for index, similar in enumerate(binary_matrix[rep_index]):
        if similar == 1:
            neighbors.append(pdbids[index])
            neighbors_index.append(index)

    new_binary_matrix = []
    new_pdbids = []
    for index1, row1 in enumerate(binary_matrix):
        row = []
        if index1 in neighbors_index: continue
        new_pdbids.append(pdbids[index1])
        for index2, similar in enumerate(row1):
            if index2 in neighbors_index: continue
            row.append(similar)
        new_binary_matrix.append(row)

    return rep, neighbors, new_binary_matrix, new_pdbids

def renumber_classes(all_classes, discrete_peps):

    new_classes = {}
    new_discrete_pep = {}

    class_count = {}

    for class_num, class_pdbids in all_classes.items():
        class_count[class_num] = len(class_pdbids)

    class_count = dict(sorted(class_count.items(), key=lambda item: item[1], reverse=True))

    new_class_num = 1
    for class_num in class_count.keys():
        new_classes[new_class_num] = all_classes[class_num]
        new_discrete_pep[new_class_num] = discrete_peps[class_num]
        new_class_num += 1

    return new_classes, new_discrete_pep

def create_greedy_classes(pdb_dict, anchor_class):

    folder = f"{DATABASE_SCRIPT_PATH}/db_stats/greedy_classification/"
    output_classes_file = f"{folder}/greedy_anchor{anchor_class}.txt"

    pdb_anchor = read_anchor_file(f"{DATABASE_SCRIPT_PATH}/db_stats/anchor_class/anchor_class.csv", anchor_class, pdb_dict)

    print(f"Greedy Classification of {len(pdb_anchor)} pHLA ∆{anchor_class}")

    pdb_dihedrals = read_central_dihedral(f"{DATABASE_SCRIPT_PATH}/db_stats/peptide_dihedrals/peptide_dihedrals_", pdb_anchor)

    binary_matrix = get_binary_matrix(pdb_dihedrals, anchor_class)

    pdbids = list(pdb_dihedrals.keys())
    all_classes = defaultdict(list)
    discrete_peps = defaultdict(str)

    class_num = 1
    while True:

        neighbors = defaultdict(list)

        rep, neighbors, binary_matrix, pdbids = get_greedy_class(binary_matrix, pdbids)
        all_classes[class_num] = neighbors
        discrete_peps[class_num] = rep
        class_num += 1

        if binary_matrix == []: break

    print(f"Classification complete (n = {len(all_classes)})")

    new_classes, new_discrete_pep = renumber_classes(all_classes, discrete_peps)

    with open(output_classes_file, "w") as txtfile:
        for class_num, class_pdbids in new_classes.items():
            txtfile.write(f"{class_num}\t{new_discrete_pep[class_num]}\t{','.join(class_pdbids)}\n")