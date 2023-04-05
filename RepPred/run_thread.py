#       Sgourakis Lab
#   Author: Sagar Gupta
#   Date: May 2, 2022
#   Email: sagarg@sas.upenn.edu

# import required libraries
import os
import argparse

import pyrosetta
pyrosetta.init()

# import custom libraries
from thread.pre_thread import PRE_THREADING
from thread.threader import THREAD
from thread.template import TEMPLATE

'''
Run threading
'''

def parse_args():

    parser = argparse.ArgumentParser(description='Adds to PDBS.fasta')

    parser.add_argument("-run_dir", help="Base directory for runs", type=str, default="threaded")
    parser.add_argument("-structure_dir", help="Directory to RELAXED structures", type=str, default=f"../")

    parser.add_argument("-target_pdbid", type=str, required=True)
    parser.add_argument("-template_pdbid", type=str, required=True)

    return parser.parse_args()

def main():

    args = parse_args()

    with open(f"./input_sequence/{args.target_pdbid}_seq.txt") as target_seqfile:
        lines = target_seqfile.readlines()
        target_hla_seq = lines[0].strip()
        target_pep_seq = lines[1].strip()

    target_seq = target_hla_seq + target_pep_seq

    pdbfile = f"{args.structure_dir}/{args.template_pdbid}_reordered.pdb"
    if not os.path.isfile(pdbfile):
        print(pdbfile, "is missing")
        quit()
    template_PDB = TEMPLATE(pdbfile)

    pre_thread = PRE_THREADING(template_PDB, target_seq, args.target_pdbid, args.run_dir)
    threader = THREAD(pre_thread)
    threader.apply()

if __name__ == "__main__":

    main()
