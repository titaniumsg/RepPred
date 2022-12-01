#       Sgourakis Lab
#   Author: Sagar Gupta
#   Date: November 29, 2022
#   Email: sagarg@sas.upenn.edu

# import required libraries
import os
from pymol import cmd

import pyrosetta
from pyrosetta import pose_from_sequence, pose_from_pdb
from pyrosetta.rosetta.core.sequence import read_aln
from pyrosetta.rosetta.protocols.comparative_modeling import PartialThreadingMover
pyrosetta.init()

'''
Threads inputted peptide and HLA sequence onto template structure
'''

def get_sequence(pdb_file, chain_list):
    
    '''
    Fetch the sequence of a PDB given the list of chain IDs
    '''
    
    seq = ''
    cmd.load(pdb_file, 'temp')
    for chain in chain_list:
        fields = cmd.get_fastastr(f"temp and chain {chain}").split("\n")
        seq += "".join(fields[1:])
    cmd.delete("all")

    return seq


def split_and_renumber(pdb_file, mhc_chain, pep_chain):
    
    '''
    Remodel spits out the PDB file with concatenated chains.

    Assumes residue 1-180 is MHC chain and residue 181- is peptide chain
    
    Renumbers residues in a structure starting from beginning
    Keeps numbering consistent with gaps if they exist in the original numbering.
    '''

    cmd.load(pdb_file)
    cmd.alter('resi 1-180', f"chain='{mhc_chain}'")
    cmd.alter("resi 181-", f"chain='{pep_chain}'")
    cmd.alter(f"chain {pep_chain}", f"resi=str(int(resi)+{-180})")
    cmd.save(pdb_file)
    cmd.delete("all")

def thread(template, output_prefix, peptide_sequence, hla_sequence):

    threaded_filename = f"{output_prefix}_{os.path.split(template)[1]}"
    return threaded_filename
    target_seq = hla_sequence+peptide_sequence
    template_seq = get_sequence(template, chain_list=['A', 'B'])  # should be 180 + 9 = 189 residues
    
    # make grishin file
    with open(f"{output_prefix}.grishin", "w") as grishinfile:
        grishinfile.write(f"## target template\n") # irrelevant
        grishinfile.write(f"#\n")
        grishinfile.write(f"scores_from_program: 0\n") # assumes no gaps in alignment
        grishinfile.write(f"0 {target_seq}\n")
        grishinfile.write(f"0 {template_seq}\n")
        grishinfile.write("--\n")

    # thread target sequence onto template structure
    template_pose = pose_from_pdb(template)
    alignment_vec = read_aln("grishin", f"{output_prefix}.grishin")
    for vec in alignment_vec:
        new_pose = pose_from_sequence(target_seq, "fa_standard")
        thread = PartialThreadingMover(vec, template_pose)
        thread.apply(new_pose)
        new_pose.dump_pdb(threaded_filename)

        split_and_renumber(threaded_filename, 'A', 'B')
    
    return threaded_filename
