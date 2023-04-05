#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: December 1, 2017
#   Email: snerli@ucsc.edu
#



# import rosetta files
import pyrosetta

from pyrosetta.rosetta.protocols.comparative_modeling import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.core.sequence import *

#custom libraries
from thread.grishin import GRISHIN


'''

PRE_THREADING class contains all the necessary functionalities required to
perform operations before the threading process.

'''

class PRE_THREADING:

    # constructor
    template = None # store the template structure
    grishin_file_name = "" # grishin file name
    target_pdbid = "" # target complex header name
    target_seq = "" # target sequence
    peptide = ""

    # constructor
    def __init__(self, template, target_seq, target_pdbid, run_dir):
        self.run_dir = run_dir
        self.template = template
        self.target_seq = target_seq
        self.target_pdbid = target_pdbid
        self.create_grishin(self.target_seq, self.template.get_sequence())

    # method to return target file name
    def get_target_file_name(self):
        return self.run_dir+"/"+self.target_pdbid+"/"+self.target_pdbid+"_on_"+self.template.get_stripped_name()

    # method to return target sequence
    def get_target_sequence(self):
        return self.target_seq

    # method to get the grishin file name
    def get_grishin_file_name(self):
        return self.grishin_file_name

    # method to return the template
    def get_template(self):
        return self.template

    # method to return the mhc header
    def get_mhc_header(self):
        return self.target_pdbid.split("_")[0]

    # method to check if grishin file exists
    def check_if_grishin_file_exists(self, filename):
        for f in self.grishin_file_name:
            if f == filename:
                return True

    # method to create a grishin file
    def create_grishin(self, target, query):
        # template_seq = query_sequence , target_seq =  target_sequence
        grishin = GRISHIN(self.get_target_file_name(), self.target_pdbid.split("_")[0],
                            self.template.get_name(), target, query)

        if not self.check_if_grishin_file_exists(grishin.get_file_name()):
            self.grishin_file_name = grishin.get_file_name()
            grishin.write()
