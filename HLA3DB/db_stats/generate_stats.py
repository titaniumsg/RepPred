#       Sgourakis Lab
#   Author: Sagar Gupta
#   Date: April 28, 2022
#   Email: sagarg@sas.upenn.edu

# import required libraries
import csv

import pyrosetta
pyrosetta.init()

from tqdm import tqdm

# import custom libraries
from misc.fasta import FASTA
from misc.constants import *
from misc.supertype import get_supertype
from misc.common import run_command, read_all_dihedrals

from db_stats.peptide_dihedrals.dihedrals import get_dihedrals
from db_stats.anchor_class.anchor_class import get_anchor_class
from db_stats.greedy_classification.greedy_classification import create_greedy_classes

import time
startTime = time.time()

def create_db():

    fasta = FASTA(f"{DATABASE_SCRIPT_PATH}/create_db/lib/PDBS.fasta")
    fasta.read()
    pdb_dict = fasta.get_pdb_dict()
    all_headers = list(fasta.get_headers())

    # get_anchor_class(pdb_dict)

    for pep_len in range(PEP_LOWER, PEP_UPPER+1):
        pdb_dict_temp = fasta.get_pdb_dict(pep_len)
        # get_dihedrals(pdb_dict_temp, pep_len)

    create_greedy_classes(pdb_dict, 6)
    create_greedy_classes(pdb_dict, 7)
    create_greedy_classes(pdb_dict, 8)

    '''

    PDB ID
    Year
    Resolution

    Peptide Sequence
    Peptide Length

    Full Allele
    Supertype
    Allele Sequence

    Anchor Class
    Anchor Distance
    Anchor1
    Anchor2

    Dihedrals

    '''
    all_dihedrals_8 = read_all_dihedrals(f"{DATABASE_SCRIPT_PATH}/db_stats/peptide_dihedrals/peptide_dihedrals_8.csv", pdb_dict)
    all_dihedrals_9 = read_all_dihedrals(f"{DATABASE_SCRIPT_PATH}/db_stats/peptide_dihedrals/peptide_dihedrals_9.csv", pdb_dict)
    all_dihedrals_10 = read_all_dihedrals(f"{DATABASE_SCRIPT_PATH}/db_stats/peptide_dihedrals/peptide_dihedrals_10.csv", pdb_dict)

    with open(f"{DATABASE_SCRIPT_PATH}/db_stats/database.csv", "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["pdbid", "year", "resolution",
                        "peptide_sequence", "peptide_length",
                        "full_allele", "supertype", "allele_sequence", "mutation",
                        "anchor_class",
                        "greedy_class",
                        "phi_4", "psi_4", "phi_5", "psi_5", "phi_6", "psi_6", "phi_7", "psi_7"])
        
        for pdbid, pep_allele in tqdm(pdb_dict.items()):

            header_info = []
            for h in all_headers:
                if h.split("|")[0] == pdbid:
                    header_info = h.split("|")
                    break

            year = int(header_info[-2])
            resolution = float(header_info[-1])

            peptide_sequence = pep_allele.split("_")[0]
            peptide_length = len(peptide_sequence)

            full_allele = pep_allele.split("_")[1]
            supertype = get_supertype(full_allele)
            if supertype == 'A_Unclassified':
                supertype = 'Unclassified A'
            elif supertype == 'B_Unclassified':
                supertype = 'Unclassified B'
            allele_sequence = run_command(f"cat {DATABASE_SCRIPT_PATH}/create_db/lib/PDBS.fasta | grep -A1 '{pdbid}|{DEFAULT_MHC_CHAIN}' | tail -1").strip()
            mutation = header_info[-3]
            if mutation == full_allele:
                mutation = 'NA'
            full_allele = f"HLA-{full_allele[0]}*{full_allele[1:3]}:{full_allele[3:]}"

            anchor_info = run_command(f"cat {DATABASE_SCRIPT_PATH}/db_stats/anchor_class/anchor_class.csv | grep '{pdbid},'").strip()
            anchor_class = int(anchor_info.split(",")[1])

            greedy_class_number = "NA"  # none should have this...
            try:
                greedy_class_info = run_command(f"cat {DATABASE_SCRIPT_PATH}/db_stats/greedy_classification/greedy_anchor{anchor_class}.txt | grep '{pdbid},'", ignore=True).strip()
                greedy_class_number = greedy_class_info.split("\t")[0]
            except:
                greedy_class_info = run_command(f"cat {DATABASE_SCRIPT_PATH}/db_stats/greedy_classification/greedy_anchor{anchor_class}.txt | grep '{pdbid}$'", ignore=True).strip()
                greedy_class_number = greedy_class_info.split("\t")[0]

            if greedy_class_number == "NA":
                print(pdbid, anchor_class)
                quit()
            
            greedy_class_number = f"{anchor_class}-{greedy_class_number}"

            dihedrals = ["NA"]*8

            if peptide_length == 8:
                dihedrals = all_dihedrals_8[pdbid][6:14]

            elif peptide_length == 9:
                dihedrals = all_dihedrals_9[pdbid][6:14]

            elif peptide_length == 10:
                dihedrals = all_dihedrals_10[pdbid][6:14]

            row = [pdbid, year, resolution,
                            peptide_sequence, peptide_length,
                            full_allele, supertype, allele_sequence, mutation,
                            anchor_class,
                   greedy_class_number] + dihedrals
            writer.writerow(row)

    executionTime = (time.time() - startTime)
    print('Script ran in', str(round(executionTime, 3)), 'seconds')

if __name__ == "__main__":

    create_db()
