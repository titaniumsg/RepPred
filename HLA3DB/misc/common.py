#       Sgourakis Lab
#   Author: Sagar Gupta
#   Date: April 28, 2022
#   Email: sagarg@sas.upenn.edu
#

# import required libraries
import subprocess
import math
import csv
from collections import defaultdict

# import custom libraries
from misc.constants import *

'''

Common methods used throughout the application

'''

def run_command(command, ignore=False):

    '''
    Helps run command line operations
    '''

    cmd = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    output = cmd.communicate()

    if cmd.returncode != 0:
        if ignore:
            return False
        print("--- FAIL ---\n")
        print(command)
        if output[0] != None:
            print(str(output[0], 'utf-8'))
        if output[1] != None:
            print(str(output[1], 'utf-8'))
        print("\n--- FAIL ---")
        return False
    else:
        if output[0] != None:
            output = output[0].decode('utf-8')
        return output
        # return True

def check_for_duplicates(query, file):

    filereader = open(file, "r")
    if query in filereader.read():
        return True
    else:
        return False

def distance(angle1, angle2):

    '''
    Determines the dihedral difference between two angles
    '''

    distance = 2.0 * (1.0 - math.cos(math.radians(angle1) - math.radians(angle2)))
    return distance

def read_all_dihedrals(dihedralsfile, pdb_dict):

    '''
    Returns phi, psi angles from the inputted dihedrals file
    '''

    pdb_dihedrals = defaultdict(list)

    with open(dihedralsfile, "r") as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for line in reader:
            pdbid = line[0]
            if pdbid not in pdb_dict.keys(): continue
            angles = line[1:]
            angles = [float(angle) for angle in angles]
            pdb_dihedrals[pdbid] = angles

    return pdb_dihedrals

def get_distance_score(pdb1, pdb2, pdb1_pepchain, pdb2_pepchain):

    

    '''
        Computes distance score between two input .PDB files 
        Only works on ∆7 9mer structures
        Assumptions:
            Both structures are 9mers
    '''

    import pyrosetta
    from pyrosetta import pose_from_pdb
    pyrosetta.init()

    pose_pdb1 = pose_from_pdb(pdb1)
    pdb1_dihedrals = get_pep_dihedrals(pose_pdb1, pdb1_pepchain, 9)
    pdb1_dihedrals_47 = pdb1_dihedrals[6:14]

    pose_pdb2 = pose_from_pdb(pdb2)
    pdb2_dihedrals = get_pep_dihedrals(pose_pdb2, pdb2_pepchain, 9)
    pdb2_dihedrals_47 = pdb2_dihedrals[6:14]

    similar, dscore = distance_score_threshold(pdb1_dihedrals_47, pdb2_dihedrals_47)

    return similar, dscore

def get_pep_dihedrals(pose_pdb, pdb1_pepchain, pep_length):
    '''
        Determines peptide dihedral angles of a given .PDB file
    '''

    dihedrals = []

    pdb1_p1 = pose_pdb.pdb_info().pdb2pose(pdb1_pepchain, 1)
    pdb1_pO = pose_pdb.pdb_info().pdb2pose(pdb1_pepchain, pep_length)

    for resi_num in range(pdb1_p1, pdb1_pO+1):
        phi = pose_pdb.phi(resi_num)
        psi = pose_pdb.psi(resi_num)

        dihedrals.append(phi)
        dihedrals.append(psi)

    return dihedrals

def distance_score_threshold(pdb1_dihedrals, pdb2_dihedrals):

    '''
        Determines if two structures are similar via distance score
        Returns (bool, float): (whether two structures are similar, their distance score)
    '''

    similar = True
    dscore = 0
    for index, angle in enumerate(pdb1_dihedrals):
        angle1 = angle
        angle2 = pdb2_dihedrals[index]

        angle_dist = distance(angle1, angle2)

        dscore += angle_dist
        if angle_dist > SINGLE_DIST_THRESHOLD:
            similar = False

    if dscore > SUM_DIST_THRESHOLD:
        similar = False

    return (similar, dscore)