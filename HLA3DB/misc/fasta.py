#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Modified: Sagar Gupta
#   Date: April 28, 2022
#   Email: snerli@ucsc.edu
#

# import other required libraries
import os
import sys
import subprocess
from collections import defaultdict

'''

FASTA class contains all the necessary functionalities required to read
and write to fasta sequence files.

'''

class FASTA:

    # class members
    filename = "" # input file name
    sequences = defaultdict(dict) # sequences read from the fasta file name

    # constructor
    def __init__(self, filename):
        self.filename = filename
        self.sequences = defaultdict(dict)

    # method to get the sequence headers
    def get_headers(self):
        if not self.sequences:
            self.sequences['none'] = ''
        return self.sequences.keys()
        
    # method to get sequences read from fasta file
    def get_sequences(self):
        return self.sequences

    # method to get a sequence given its header
    def get_sequence(self, header):
        return self.sequences[header]

    # method to read the fasta file
    # populates headers and sequences
    def read(self):
        if self.filename != None:
            fasta_file_handler = open(self.filename, "r")
            try:
                for line in fasta_file_handler:
                    line = line.rstrip()
                    if ">" in line:
                        header = "".join(list(line)[1:])
                        self.sequences[header] = ''
                    else:
                        self.sequences[header] += line
                fasta_file_handler.close()
            except UnicodeDecodeError:
                print ("UnicodeDecodeError")
                fasta_file_handler.close()

    # method to get dictionary
    # of the format {PDB_ID: "PeptideSequence_Allele"}
    def get_pdb_dict(self, pep_len=None):
        if self.filename != None:
            fasta_file_handler = open(self.filename, "r")
            counter = 1
            raw_dict = {}
            pdbid = ""
            for line in fasta_file_handler:
                if ">" in line:
                    line = line.rstrip().lstrip('>')
                    pdbid = line.split('|')[0]
                    allele = line.split('|')[2]
                    if counter % 2 == 0:
                        raw_dict[pdbid] = ""
                    counter += 1
                else:
                    line = line.rstrip()
                    if pep_len != None:
                        if len(list(line)) == pep_len:
                            raw_dict[pdbid] = line+"_"+allele
                    else:
                        raw_dict[pdbid] = line+"_"+allele
            fasta_file_handler.close()
            pdb_dict = {}
            for key, value in raw_dict.items():
                if value != "":
                    pdb_dict[key] = value

            return pdb_dict

    def get_pdb_dict_year(self, pep_len=None):
        if self.filename != None:
            fasta_file_handler = open(self.filename, "r")
            pdb_dict_year = {}
            for line in fasta_file_handler:
                if ">" in line:
                    line = line.rstrip().lstrip('>')
                    pdbid = line.split('|')[0]
                    year = int(line.split('|')[-2])
                    pdb_dict_year[pdbid] = year
            fasta_file_handler.close()

            return pdb_dict_year

    # method to write to fasta file
    def write(self, header, text):
        fasta_file_handler = open(self.filename, "w")
        fasta_file_handler.write(">"+header+"\n")
        fasta_file_handler.write(text+"\n")
        fasta_file_handler.close()

    def get_pdb_dict_allele(self, selected_allele, pep_len=None):
        if self.filename != None:
            fasta_file_handler = open(self.filename, "r")
            counter = 1
            raw_dict = {}
            pdbid = ""
            for line in fasta_file_handler:
                if ">" in line:
                    line = line.rstrip().lstrip('>')
                    pdbid = line.split('|')[0]
                    allele = line.split('|')[2]
                    if counter % 2 == 0:
                        raw_dict[pdbid] = ""
                    counter += 1
                else:
                    line = line.rstrip()
                    if pep_len != None:
                        if len(list(line)) == pep_len:
                            raw_dict[pdbid] = line+"_"+allele
                    else:
                        raw_dict[pdbid] = line+"_"+allele
            fasta_file_handler.close()
            pdb_dict = {}
            for key, value in raw_dict.items():
                if value != "" and value.split("_")[1] == selected_allele:
                    pdb_dict[key] = value

            return pdb_dict
