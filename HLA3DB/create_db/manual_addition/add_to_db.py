#       Sgourakis Lab
#   Author: Sagar Gupta
#   Date: April 28, 2022
#   Email: sagarg@sas.upenn.edu

# import required libraries
import glob
import os

import Bio
from Bio.PDB import PDBParser

# import custom libraries
from misc.constants import *
from misc.common import run_command

from create_db.update_db import get_pdb, get_sequence, run_align, get_hla_sequences, is_seq_present, rename_chains, fetch_dates_of_pdbs, update_PDBS_fasta

def add_to_db():

    for path_pdbfile in glob.glob(f"{DATABASE_SCRIPT_PATH}/create_db/manual_addition/structures/*.pdb"):
        pdbfile = os.path.basename(path_pdbfile)
        pdbid = pdbfile.split("_reordered")[0]

        raw_pdbid = pdbid[0:4]

        _, resolution = get_pdb(raw_pdbid)
        run_command("rm -r temp1 obsolete")

        mhc_chains, peptide_chains = [], []
        mhc_seqs = []

        p = PDBParser()
        s = p.get_structure('X', path_pdbfile)

        for model in s:
            for ch in model:
                chain = ch.id
                seq = get_sequence(path_pdbfile, chain)

                if len(seq) >= PEP_LOWER and len(seq) <= PEP_UPPER:
                    peptide_chains.append(chain)
                    pep_length = len(seq)

                else:
                    (align1, align2, no_gaps, score) = run_align(A0201_SEQ, seq)
                    if score > HC_CUTOFF:
                        mhc_chains.append(chain)
                        mhc_seqs.append(seq)
                    else:
                        continue

        final_mhc_chain = mhc_chains[0]
        final_pep_chain = peptide_chains[0]
        final_mhc_seq = mhc_seqs[0]

        seq2allele = get_hla_sequences()
        (first_allele, raw_allele) = is_seq_present(seq2allele, get_sequence(path_pdbfile, final_mhc_chain))

        mhc_chain, pep_chain = rename_chains(path_pdbfile, final_mhc_chain, final_pep_chain)

        allele_seq = get_sequence(path_pdbfile, mhc_chain)
        pep_seq = get_sequence(path_pdbfile, pep_chain)
        pMHC = pep_seq+"_"+first_allele

        pdbid_date = fetch_dates_of_pdbs(jsonfile=f"{DATABASE_SCRIPT_PATH}/create_db/rcsb_search.json", debug_pdbid=None)
        if raw_pdbid in pdbid_date.keys():
            pdbid_date[pdbid] = pdbid_date[raw_pdbid]

        update_PDBS_fasta(pdbid, raw_allele, mhc_chain, pep_chain, allele_seq, pep_seq, f"{DATABASE_SCRIPT_PATH}/create_db/lib/PDBS.fasta", pdbid_date, resolution)
        run_command(f"cp -v {path_pdbfile} {DATABASE_SCRIPT_PATH}/create_db/lib/TRIMMED_TEMPLATES/")
    
    run_command("rm rcsb_search_temp.json")


if __name__ == "__main__":

    add_to_db()
