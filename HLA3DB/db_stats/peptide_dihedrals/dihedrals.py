#       Sgourakis Lab
#   Author: Sagar Gupta
#   Date: April 28, 2022
#   Email: sagarg@sas.upenn.edu

# import required libraries
import csv

import pyrosetta
from pyrosetta import pose_from_pdb

# import custom libraries
from misc.constants import *
from misc.common import check_for_duplicates

def get_dihedrals(pdb_dict, pep_len):

    outputfile = f"{DATABASE_SCRIPT_PATH}/db_stats/peptide_dihedrals/peptide_dihedrals_{pep_len}.csv"

    with open(outputfile, "a") as dihedralsfile:
        writer = csv.writer(dihedralsfile)

        header = ["pdbid"]
        for i in range(1, pep_len+1):
            header.append(f"phi_{i}")
            header.append(f"psi_{i}")

        if not check_for_duplicates(",".join(header), outputfile):
            writer.writerow(header)

        for pdbid in pdb_dict.keys():
            row = [pdbid]
            pose = pose_from_pdb(f"{DATABASE_SCRIPT_PATH}/create_db/lib/TRIMMED_TEMPLATES/{pdbid}_reordered.pdb")
            p1 = pose.pdb_info().pdb2pose(DEFAULT_PEP_CHAIN, 1)
            pO = pose.pdb_info().pdb2pose(DEFAULT_PEP_CHAIN, pep_len)

            for resi_num in range(p1,pO+1):
                phi = pose.phi(resi_num)
                psi = pose.psi(resi_num)

                row.append(str(phi))
                row.append(str(psi))

            if not check_for_duplicates(",".join(row), outputfile):
                writer.writerow(row)
