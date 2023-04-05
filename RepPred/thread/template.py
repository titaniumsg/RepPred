#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Date: January 26, 2017
#   Email: snerli@ucsc.edu

# Load the Rosetta commands for use in the Python shell
from pyrosetta import *

# import rosetta files
from pyrosetta.rosetta.protocols.relax import *
from pyrosetta.rosetta.core.scoring import *

'''

TEMPLATE class contains all the necessary functionalities required to
perform prepare the template structure for threading.

'''


class TEMPLATE:

    # class members
    template_pdb = None # store template pdb file name
    template_pose = None # store the template pose

    # constructor
    def __init__(self, pdb):
        self.template_pdb = pdb

    # method to return pose
    def get_pose(self):
        return self.template_pose

    # method to return sequence of the template
    def get_sequence(self):
        if self.template_pose == None:
            self.get_pose_from_pdb()
        return self.template_pose.sequence()

    # method to return name of the template pdb file
    def get_pdb(self):
        return self.template_pdb

    # method to load pose from a pdb file
    def get_pose_from_pdb(self):
        if self.template_pose == None:
            self.template_pose = pose_from_pdb(self.template_pdb)
        return self.template_pose

    # method to return pdb name fields joined with dot
    # some of the naming convention we use
    def get_name(self):
        template_pdb_fields = self.template_pdb.split(".")
        template_pdb_name = ".".join(template_pdb_fields[0:len(template_pdb_fields)-1])
        return template_pdb_name

    # method to get the striped pdb name or file name without the .pdb extension
    def get_stripped_name(self):
        template_pdb_name = os.path.basename(self.template_pdb)
        return template_pdb_name.split('.')[0]

    # methot to return absolute path of the file where template is present
    def get_template_path(self):
        return os.path.abspath(os.path.dirname(self.template_pdb))+"/"

    # method to save pdb file
    def save_pdb(self):
        self.template_pose.dump_pdb(self.template_pdb)
