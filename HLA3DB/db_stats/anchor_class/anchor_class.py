#       Sgourakis Lab
#   Author: Sagar Gupta
#   Date: April 28, 2022
#   Email: sagarg@sas.upenn.edu

# import required libraries
import argparse
import csv

import pyrosetta
from pyrosetta import pose_from_pdb

# import custom libraries
from misc.constants import *

def get_ca_ca_distance(pose, resi1, resi2):
    return pose.residue(resi1).atom("CA").xyz().distance(pose.residue(resi2).atom("CA").xyz())

def get_anchor_class(pdb_dict):

    with open(f"{DATABASE_SCRIPT_PATH}/db_stats/anchor_class/anchor_class.csv", "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["pdbid", "anchor_class", "anchor_dist", "anchor1", "anchor2"])
        for pdbid, pep_allele in pdb_dict.items():
            pep_len = len(pep_allele.split("_")[0])

            pose = pose_from_pdb(f"{DATABASE_SCRIPT_PATH}/create_db/lib/TRIMMED_TEMPLATES/{pdbid}_reordered.pdb")
            p1 = pose.pdb_info().pdb2pose(DEFAULT_PEP_CHAIN, 1)
            pO = pose.pdb_info().pdb2pose(DEFAULT_PEP_CHAIN, pep_len)

            dist_1 = {}
            dist_2 = {}

            for resi_num in range(p1,pO+1):
                dist_1[resi_num] = get_ca_ca_distance(pose, 24, resi_num)
                dist_2[resi_num] = get_ca_ca_distance(pose, 123, resi_num)

            anchor_1 = min(dist_1, key=dist_1.get)
            anchor_2 = min(dist_2, key=dist_2.get)
            anchor_class = anchor_2 - anchor_1

            anchor_dist = get_ca_ca_distance(pose, anchor_1, anchor_2)

            writer.writerow([pdbid, anchor_class, anchor_dist, anchor_1, anchor_2])
