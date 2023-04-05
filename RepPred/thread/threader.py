#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: December 1, 2017
#   Email: snerli@ucsc.edu
#

# Load the Rosetta commands for use in the Python shell

# import rosetta files
import pyrosetta
from pyrosetta import pose_from_sequence

from pyrosetta.rosetta.protocols.comparative_modeling import PartialThreadingMover
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.core.sequence import *

#custom libraries

# import other required libraries
import pymol
from Bio.PDB import PDBParser, PDBIO
'''

THREAD class contains all the necessary functionalities required to
perform target onto the template structure.

'''

class THREAD:

    # class members
    pre_threader = None # pre threader object
    post_threader = None
    threaded_pose = None
    tag = None

    # constructor
    def __init__(self, pre_thread):
        self.pre_threader = pre_thread

    # method that calls threader
    def apply(self):
        if self.pre_threader != None:
            self.do_partial_threading()

    # method to perform partial threading
    def do_partial_threading(self):
        # for each alignment combination in the grishin file.
        alignment_vec = read_aln("grishin", self.pre_threader.get_grishin_file_name())
        for vec in alignment_vec:
            new_pose = pose_from_sequence(self.pre_threader.get_target_sequence(), "fa_standard")
            thread = PartialThreadingMover(vec, self.pre_threader.get_template().get_pose())

            thread.apply(new_pose)

            self.threaded_pose = new_pose
            self.tag = self.pre_threader.get_target_file_name()

            # output threaded structure
            self.threaded_pose.dump_pdb(self.tag+".pdb")
            self.split_pep_mhc_chains()
            self.renumber_residues()

    def split_pep_mhc_chains(self):

        target_file = self.tag+".pdb"
        pymol.cmd.load(target_file)

        pymol.cmd.do("select resi 181-")
        pymol.cmd.extract("pep", "sele")
        pymol.cmd.alter('pep and chain A', 'chain="B"')
        pymol.cmd.create('merged', "all")
        pymol.cmd.do("save "+self.tag+".pdb"+", merged")
        pymol.cmd.do("delete all")

    def renumber_residues(self, begin=1):
        """ Renumbers residues in a structure starting from begin (default: 1).
            Keeps numbering consistent with gaps if they exist in the original numbering.
        """

        p = PDBParser()
        structure = p.get_structure('X', self.tag+'.pdb')
        for model in structure:
            for chain in model:
                if chain.get_list()[0].get_id()[1] == 0:
                    return structure
                fresidue_num = chain.get_list()[0].get_id()[1]
                displace = begin - fresidue_num
                for res in chain:
                    res.id = (res.id[0], res.id[1]+displace, res.id[2])

        io = PDBIO()
        io.set_structure(structure)
        io.save(self.tag+'.pdb')
