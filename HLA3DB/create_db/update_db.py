#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Modified: Sagar Gupta
#   Date: August 12, 2022
#   Email: snerli@ucsc.edu and sagarg@upenn.edu

# import required libraries
import os
import json
import math
import time
import requests
import argparse
import warnings
import platform
from collections import defaultdict
from datetime import date

import Bio
from Bio.PDB import *
from Bio.PDB import PDBList
from Bio.PDB import PDBParser
from Bio.PDB import Select, PDBIO
from Bio import pairwise2

import pymol
from pymol import cmd, stored

# import custom libraries
from misc.HLA_sequences_180 import hla_sequences_180
from misc.fasta import FASTA
from misc.constants import *
from misc.common import run_command

# As of November 2022, script takes ~66 minutes to run
startTime = time.time()
warnings.filterwarnings("ignore")


'''
Script to auto-update the database in HLA3DB

python update_db.py -rosettainstalldir <ROSETTA_FOLDER>
'''


def get_args():
    """
    Parse command line arguments
    """

    parser = argparse.ArgumentParser(description="Method to update database with latest pHLA-I structures in the PDB")

    parser.add_argument("-rosettainstalldir", help="Rosetta install directory to run RosettaRemodel for missing N terminal residues", default=ROSETTA_INSTALL_DIR)

    parser.add_argument("-start_fasta", help="fasta file of starting list of structures, make file blank if starting from scratch", default="PDBS_empty.fasta")
    parser.add_argument("-date", help="If date (year) is not specified, the script will start from 1988", default=None)

    parser.add_argument("-debug_pdbid", help="PDB ID for debugging purpose only")

    args = parser.parse_args()

    return args

def get_hla_sequences():
    """
    Get HLA sequences from the database in HLA_sequences_180.py
    and return the reversed key-value pairs of the hla_sequences_180 dictionary
    """

    seq2allele = defaultdict(list)
    for key, value in hla_sequences_180.items():
        seq2allele[value].append(key)
    return seq2allele


def update_json(jsonfile, date_to_fetch):
    """
    Update and create JSON for querying PDB based on release date of the structure
    Simply replaces the year (default 1988) with the inputted year (date_to_fetch)
    """

    if date_to_fetch != None:
        outputfilehandler = open(f"{jsonfile.split('.json')[0]}_temp.json", "w")
        inputfilehandler = open(jsonfile, 'r')
        for line in inputfilehandler:
            if "1988" in line:
                outputfilehandler.write(line.replace("1988", date_to_fetch))
            else:
                outputfilehandler.write(line)
        inputfilehandler.close()
        outputfilehandler.close()


def fetch_dates_of_pdbs(jsonfile, debug_pdbid=None):
    """
    Fetch release dates of PDBIDs
    """

    pdbid_date = {}
    if debug_pdbid == None:

        # assumes you want it from the year you are in
        todays_year = date.today().year

        rcsb = 'https://search.rcsb.org/rcsbsearch/v2/query'
        for year in range(FIRST_MHC_STRUC_YEAR, todays_year+2):

            update_json(jsonfile, str(year))
            jsonfilehandler = open(f"{jsonfile.split('.json')[0]}_temp.json",)
            payload = json.load(jsonfilehandler)
            rcsb_output = requests.post(rcsb, json=payload)
            if rcsb_output.status_code != 200:
                continue
            rcsb_data = rcsb_output.json()
            num_structures = int(rcsb_data['total_count'])

            for i in range(0, num_structures):
                pdbid_date[rcsb_data['result_set'][i]['identifier']] = str(year)
    else:
        pdbid_date[debug_pdbid] = "1950"  # set debug year as 1950 for easy identification
    return pdbid_date


def fetch_pdbids(jsonfile, date, debug_pdbid=None):
    """
    Fetch PDBIDs of the structures deposited since a specified date
    """

    pdbids = []

    if debug_pdbid == None:

        rcsb = 'https://search.rcsb.org/rcsbsearch/v2/query'
        if date == None:
            jsonfilehandler = open(jsonfile,)
        else:
            jsonfilehandler = open(f"{jsonfile.split('.json')[0]}_temp.json",)
        payload = json.load(jsonfilehandler)
        rcsb_output = requests.post(rcsb, json=payload)
        if rcsb_output.status_code != 200:
            return pdbids
        rcsb_data = rcsb_output.json()
        num_structures = int(rcsb_data['total_count'])
        for i in range(0, num_structures):
            pdbids.append(rcsb_data['result_set'][i]['identifier'])

    else:
        pdbids.append(debug_pdbid)

    return pdbids


def get_unique_pdbids(fastafile, pdbids):
    """
    Check if the PDBIDs are unique relative to the earlier manually curated PDB list
    """

    fasta = FASTA(fastafile)
    fasta.read()

    headers = fasta.get_headers()
    pdbs_curated = []
    pep_seq_allele_pair = {}
    for header in headers:
        fields = header.split('|')
        if len(fields) >= 2:
            pdbid = fields[0]
            allele = fields[2]
            pdbs_curated.append(pdbid)
            seq = fasta.get_sequence(header)
            if len(seq) >= PEP_LOWER and len(seq) <= PEP_UPPER:
                pep_seq_allele_pair[seq+'_'+allele] = pdbid

    pdbs_to_curate = []
    for pdbid in pdbids:
        if pdbid not in pdbs_curated:
            pdbs_to_curate.append(pdbid)

    return (pdbs_to_curate, pep_seq_allele_pair)


def get_pdb(pdbid):
    """
    Fetch the PDB and its resolution from the server given its PDB ID
    """

    try:
        pdblist = PDBList()

        pdbfile = "temp1/"+pdbid.lower()+".pdb"
        if not os.path.exists(pdbfile):
            pdblist.retrieve_pdb_file(pdbid, pdir="temp0", file_format="pdb")
            outfile = "temp0/pdb"+pdbid.lower()+".ent"
            if not os.path.exists("temp1"):
                run_command("mkdir temp1")

            run_command(f"cp {outfile} {pdbfile}")
            parser = PDBParser()
            structure = parser.get_structure(pdbid, f"temp1/{pdbid.lower()}.pdb")
            resolution = str(structure.header["resolution"])

            run_command("grep ^ATOM "+outfile+" > "+pdbfile)
            run_command("rm -r temp0")

        return pdbfile, resolution
    except:
        print(pdbid+" was not found or had some issues while downloading")
        return None, None


def get_sequence(targetfile, chainid):
    """
    Fetch the sequence of a PDB given the chain ID
    """

    targetfile_name = os.path.basename(targetfile).split(".pdb")[0]
    cmd.load(targetfile, targetfile_name)
    raw_seq = pymol.cmd.get_fastastr(f"{targetfile_name} and chain {chainid}")
    cmd.delete("all")

    fields = raw_seq.split("\n")
    seq = "".join(fields[1:])

    return seq


def run_align(refseq, targetseq, find_mutations=False):
    '''
    Aligns targetseq to refseq and outputs either the alignment and score
    or the mutations and score
    '''

    alignments = pairwise2.align.globalms(refseq, targetseq, 1, 0, -1, 0)
    no_gaps = 1000
    score = -1
    align1 = ''
    align2 = ''

    for alignment in alignments:
        if float(alignment[2]) > score:
            align1 = alignment[0]
            align2 = alignment[1]
            no_gaps = alignment[1].count('-')
            score = float(alignment[2])

    if not find_mutations:
        return (align1, align2, no_gaps, score)

    else:
        mutations = []
        for i in range(0, len(align1)):
            if align1[i] != '-' and align2[i] != '-' and align1[i] != align2[i]:
                mutations.append(align2[i]+str(i+1)+align1[i])

        return (mutations, score)


def find_max_index(arr):
    '''
    Finds maximum value and respective index of inputted list
    '''

    max_val = arr[0]
    max_index = 0
    for i in range(1, len(arr)):
        if int(max_val) < int(arr[i]):
            max_index = i
            max_val = arr[i]

    return (max_val, max_index)


def get_start_resid_for_chains(pdbfile):
    '''
    Returns a dictionary of the first residue number for each chain (usually 1)
    '''

    p = PDBParser()
    s = p.get_structure('X', pdbfile)
    resid_for_chains = {}

    for model in s:
        for ch in model:
            for r in ch:
                resid_for_chains[ch.id] = r.id[1]
                break

    return resid_for_chains


def align_chains_to_reference(pdbfile):
    '''
    Determines identity of each chain using sequence alignment to HLA-A*02:01 MHC and B2M (see constants.py)
    Also determines if structure needs to be remodeled because of incorrect HLA heavy chain N-terminal residues
    '''

    p = PDBParser()
    s = p.get_structure('X', pdbfile)
    peptide_chains = []
    mhc_chains = []
    alignment_scores = []
    b2m_chains = []
    other_chains = []
    resid_for_chains = {}

    mhc_seqs = []
    pep_length = 0
    pep_over_10AA = False

    pdb_not_hla = False
    pdb_has_receptor = False
    remodel = False
    residues_to_insert = []
    mhc_seq_truncated = False

    for model in s:
        for ch in model:
            chain = ch.id
            seq = get_sequence(pdbfile, chain)

            # determine peptide chain by length of chain
            if len(seq) >= PEP_LOWER and len(seq) <= PEP_UPPER:
                peptide_chains.append(chain)
                pep_length = len(seq)

            # determine long peptide as between 10 and 20
            elif len(seq) > PEP_UPPER and len(seq) <= LONG_PEP_LENGTH:
                pep_over_10AA = True
            else:
                (align1, align2, no_gaps, score) = run_align(A0201_SEQ, seq)

                # first check if heavy chain is valid
                if score > HC_CUTOFF:
                    mhc_chains.append(chain)
                    alignment_scores.append(score)
                    mhc_seqs.append(seq)

                # if it's not valid, check if b2m is valid
                else:
                    (align1, align2, no_gaps, score) = run_align(B2M_SEQ, seq)
                    if score > B2M_CUTOFF:
                        b2m_chains.append(chain)
                    else:
                        other_chains.append(chain)

    resid_for_chains = get_start_resid_for_chains(pdbfile)
    final_mhc_chain = ''
    final_pep_chain = ''
    ignore_receptors = False

    if len(mhc_chains) > 0 and len(peptide_chains) > 0:
        (max_score, max_index) = find_max_index(alignment_scores)

        # chooses mhc chain with the highest alignment score to A0201
        final_mhc_chain = mhc_chains[max_index]
        final_pep_chain = peptide_chains[0]
        found_pep_chain = False

        # find the best MHC chain's respective peptide chain
        for model in s:
            for p in peptide_chains:

                # 65th residue of MHC Ca atom and 2nd residue for peptide Ca atom (from observation)
                distance = model[final_mhc_chain][resid_for_chains[final_mhc_chain]+64]['CA'] - model[p][resid_for_chains[p]+1]['CA']
                if distance <= DISTANCE_PEP_MHC:
                    final_pep_chain = p
                    found_pep_chain = True
                    break
            if found_pep_chain:
                break

        if len(b2m_chains) > 0:
            if len(other_chains) > 0:
                pdb_has_receptor = True

                if not is_receptor_near_peptide(pdbfile, mhc_chains, b2m_chains, peptide_chains, other_chains, final_pep_chain):
                    ignore_receptors = True

        if not pdb_not_hla and ((not pdb_has_receptor) or (pdb_has_receptor and ignore_receptors)) and len(final_mhc_chain) > 0 and len(final_pep_chain) > 0:
            (remodel, residues_to_insert) = identify_if_missing_N_terminal_residues(pdbfile, final_mhc_chain, mhc_seqs[max_index])

        if (MHC_TRIM_LENGTH - 1) - len(residues_to_insert) - len(mhc_seqs[max_index]) > 0:
            mhc_seq_truncated = True

    elif len(other_chains) > 0:
        pdb_not_hla = True

    return (pdb_has_receptor, pdb_not_hla, final_mhc_chain, final_pep_chain, pep_length, remodel, residues_to_insert, pep_over_10AA, b2m_chains, ignore_receptors, mhc_seq_truncated)


def is_receptor_near_peptide(pdbfile, mhc_chains, b2m_chains, pep_chains, other_chains, final_pep_chain):
    '''
    Determines if a receptor (TCR, etc.) is within user-defined (default = 5 angstroms)
    distance from peptide chain
    '''

    if len(mhc_chains) > 0 and len(b2m_chains) > 0 and len(pep_chains) > 0:

        o_chains = '+'.join(other_chains)

        cmd.load(pdbfile)

        cmd.select(f"chain {final_pep_chain}")
        cmd.do("set_name sele, pep")
        cmd.select(f"neighbor_res", f"(br. pep around {str(RECEPTOR_DIST)}) and chain {o_chains}")

        stored.list = []
        cmd.iterate("(neighbor_res and name ca)", "stored.list.append((resi))")

        cmd.delete("all")

        if len(stored.list) > 0:
            return True

    else:
        return True

    return False


def identify_if_missing_N_terminal_residues(pdbfile, mhc_chain, mhc_seq):
    '''
    Determines which MHC heavy chain N-terminal residues are missing
    on the basis of the SHS motif
    '''

    residues_to_insert = []
    # first 50 residues of the aligned portion of MHC seq
    first_res = ''.join(mhc_seq[:50])

    remodel = False
    if MHC_SEQ_MOTIF_SHS in first_res:
        index = first_res.index(MHC_SEQ_MOTIF_SHS)
        # need to add G if seq begins with SHS
        if index == 0:
            remodel = True
            residues_to_insert.append('G')
        # if it is later on then we need to trim up til the GSHS part (index - 1)
        elif index > 1:
            length_to_trim = index - 1
            trim_left_end(pdbfile, mhc_chain, length_to_trim)

    elif MHC_SEQ_MOTIF_HS in first_res:
        index = first_res.index(MHC_SEQ_MOTIF_HS)
        # need to add GS if seq begins with HS
        if index == 0:
            remodel = True
            residues_to_insert.append('G')
            residues_to_insert.append('S')

    return (remodel, residues_to_insert)


def trim_left_end(pdbfile, chain, length_to_trim):
    """
    Trim the N terminal end of PDB which consists of the MHC heavy chain based on the length specified
    """

    p = PDBParser()
    s = p.get_structure('X', pdbfile)

    resid_for_chains = get_start_resid_for_chains(pdbfile)
    chain_length = resid_for_chains[chain] + length_to_trim

    class ResSelect(Select):
        def accept_residue(self, res):
            if res.parent.id == chain and res.id[1] < chain_length:
                return False
            else:
                return True

    io = PDBIO()
    io.set_structure(s)
    outfile = pdbfile.split('.pdb')[0]+"_N_trim.pdb"
    io.save(outfile, ResSelect())

    run_command("mv "+outfile+" "+pdbfile)


def run_remodel(pdbfile, mhc_chain, pep_chain, residues_to_insert, rosettainstalldir):
    """
    Remodel the N-terminal region of the heavy chain in the given PDB
    """

    if len(residues_to_insert) > 0:
        outputfilehandler = open("test.remodel", "w")
        for res in residues_to_insert:
            outputfilehandler.write(f"1 X E PIKAA {res}\n")
        outputfilehandler.close()

    # this command may need to be changed based on Rosetta version
    run_command(f"{rosettainstalldir}/main/tools/remodel/getBluePrintFromCoords.pl -pdbfile {pdbfile} -chain {mhc_chain} >> test.remodel")

    outputfilehandler = open("test_updated.remodel", "w")
    inputfilehandler = open("test.remodel", "r")
    not_start = False
    for line in inputfilehandler:
        if len(residues_to_insert) > 0:
            if "PIKAA" not in line and not not_start:
                line = line.rstrip()
                fields = line.split()
                id = fields[0]
                aa = fields[1]
                newline = id+" "+aa+" E PIKAA "+aa+"\n"
                outputfilehandler.write(newline)
                not_start = True
            else:
                outputfilehandler.write(line)
        else:
            outputfilehandler.write(line)

    inputfilehandler.close()
    outputfilehandler.close()

    os_system = platform.system()
    extension = ''
    if os_system == 'Linux':
        extension = 'linuxgccrelease'
    elif os_system == 'Darwin':
        extension = 'macosclangrelease'

    # this command may need to be changed based on Rosetta version
    # remodel.static instead of remodel depending on user files
    run_command(f"{rosettainstalldir}/main/source/bin/remodel.{extension} -database {rosettainstalldir}/main/database -s {pdbfile} -remodel:blueprint test_updated.remodel -remodel::num_trajectory 1")
    # run_command(f"{rosettainstalldir}/main/source/bin/remodel.mpi.{extension} -database {rosettainstalldir}/main/database -s {pdbfile} -remodel:blueprint test_updated.remodel -remodel::num_trajectory 1")

    new_pdbfile = pdbfile.replace(".pdb", "_0001.pdb")
    filename_fields = new_pdbfile.split('/')
    new_pdbfile_name = filename_fields[len(filename_fields)-1]

    break_and_renumber(new_pdbfile_name, len(residues_to_insert), mhc_chain, pep_chain)

    run_command(f"cp {new_pdbfile_name} {pdbfile}")


def break_and_renumber(pdbfile, gaps, mhc_chain, pep_chain):
    """
    Remodel sometimes spits out the PDB file with concatenated chains.
    We want to break the chains and renumber them
    """

    p = PDBParser()
    s = p.get_structure('X', pdbfile)

    class ResSelect(Select):
        def accept_residue(self, res):
            if res.id[1] < MHC_TRIM_LENGTH or res.id[1] >= MHC_TRIM_LENGTH:
                return True
            else:
                return False
    io = PDBIO()
    io.set_structure(s)
    outfile = pdbfile.split('.pdb')[0]+"_C_trim.pdb"
    io.save(outfile, ResSelect())

    split_pep_mhc_chains(outfile, mhc_chain, pep_chain)
    renumber_residues(outfile)

    run_command(f"cp {outfile} {pdbfile}")


def split_pep_mhc_chains(targetfile, mhc_chain, pep_chain):
    """
    Remodel spits out the PDB file with concatenated chains.
    Assumes residue 1-180 is MHC chain and residue 181- is peptide chain
    We want to break the chains
    """

    cmd.load(targetfile)
    cmd.alter('resi 1-180', f"chain='{mhc_chain}'")
    cmd.alter("resi 181-", f"chain='{pep_chain}'")
    cmd.save(targetfile)
    cmd.delete("all")


def renumber_residues(targetfile, begin=1):
    """
    Renumbers residues in a structure starting from begin (default: 1).
    Keeps numbering consistent with gaps if they exist in the original numbering.
    """

    start_resid = get_start_resid_for_chains(targetfile)
    cmd.load(targetfile)
    chains = cmd.get_chains()

    for chain in chains:
        start = begin - start_resid[chain]
        cmd.alter(f"chain {chain}", f"resi=str(int(resi)+{start})")

    cmd.save(targetfile)
    cmd.delete("all")


def is_missing_bb_heavy_atoms(filename, target_chain):
    """
    Check if a PDB is missing any backbone heavy atoms
    """

    countN = 0
    countCA = 0
    countC = 0
    countO = 0
    no_residues = 0
    missing_residue = False
    missing_heavy_atoms = False

    parser = PDBParser()
    structure = parser.get_structure('X', filename)
    for model in structure:
        for chain in model:
            if chain.id == target_chain:
                no_residues += len(get_sequence(filename, chain.id))
                prev_res_index = 0
                ind_set = False
                offset = 0
                for residue in chain:
                    if not ind_set:
                        prev_res_index = residue.id[1] - 1
                        ind_set = True
                        offset = residue.id[1] - 1
                    for atom in residue:
                        atom_name = atom.get_name()
                        if atom_name == 'N':
                            countN += 1
                        elif atom_name == 'C':
                            countC += 1
                        elif atom_name == 'CA':
                            countCA += 1
                        elif atom_name == 'O':
                            countO += 1

                    if residue.id[1] == prev_res_index+1:
                        prev_res_index = residue.id[1]
                    else:
                        missing_residue = True

                if (prev_res_index-offset) != len(get_sequence(filename, chain.id)):
                    missing_residue = True

    if countN != no_residues or countC != no_residues or countCA != no_residues or countO != no_residues:
        missing_heavy_atoms = True

    return (missing_residue or missing_heavy_atoms)


def has_zero_occupancy_atoms(pdbfile, target_chain):
    '''
    Check if there are any zero occupancy backbone heavy atoms in residues of the target chain
    '''

    parser = PDBParser()
    structure = parser.get_structure('X', pdbfile)
    for model in structure:
        for chain in model:
            # only look at peptide chain
            if chain.get_id() == target_chain:
                for residue in chain:
                    for atom in residue:
                        atom_name = atom.get_name()
                        occupancy = atom.get_occupancy()
                        # only look at peptide backbone heavy atoms
                        if atom_name == 'N' or atom_name == 'C' or atom_name == 'CA' or atom_name == 'O':
                            # zero occupancy means the location and properties of the residue may not be reliable
                            if occupancy == 0.0:
                                print(pdbfile, 'has an unreliable peptide backbone at', str(residue.id[1])+atom_name)
                                return True

    return False


def is_seq_present(seq2allele, seq):
    """
    Check if a sequence is present in the allele database
    """

    for key, val in seq2allele.items():
        if str(seq) in key:
            raw_allele = sorted(val)[0]
            return (raw_allele.replace("*", "").replace(":", ""), raw_allele)

    score = -1
    mutations = []
    allele = ''
    raw_allele = ''
    for key, val in seq2allele.items():
        (mutations, new_score) = run_align(key, seq, True)
        if int(new_score) > int(score) and len(mutations) <= MUTATION_THRESH:
            raw_allele = sorted(val)[0]+'|'+','.join(mutations)
            allele = raw_allele.replace("*", "").replace(":", "")+'|'+','.join(mutations)
            score = new_score

    return (allele, raw_allele)


def distance(angle1, angle2):
    """
    Function to determine the difference between two dihedral angles (phi or psi)
    """

    distance = 2.0 * (1.0 - math.cos(math.radians(angle1) - math.radians(angle2)))
    return distance


def resolution_of_bb(pdbfile, pep_chain, pep_length):
    """
    Checks if there is only one peptide backbone within a certain threshold (18.2 degrees)
    """

    cmd.load(pdbfile)
    phi_psi = cmd.phi_psi("chain "+pep_chain)
    cmd.delete("all")

    # quick check to see if number of angles is exact (e.g. 9mer has 14 dihedrals, if above 14 -> checks for double backbone)
    if len(phi_psi) == (int(pep_length) - 2):
        return True
    else:
        resi_num = []
        angles = []
        for pdb, dihedrals in phi_psi.items():
            resi_num.append(pdb[1])
            angles.append(dihedrals)

        for i in range(0, len(resi_num)):
            # checks neighboring "residue"
            if (resi_num[i] + 1) in resi_num:
                # checks distance between phi angles
                low_res_bb_dist_phi = distance(angles[i][0], angles[i+1][0])
                # checks distance between psi angles
                low_res_bb_dist_psi = distance(angles[i][1], angles[i+1][1])
                # checks if both phi and psi distances are under a threshold (18.2 degrees)
                if low_res_bb_dist_phi < 0.1 and low_res_bb_dist_psi < 0.1:
                    continue
                else:
                    return False
    return True


def trim_and_save_mhc_pep_chains(pdbfile, pdbid, mhc_chain, pep_chain, pep_length, residues_to_insert):
    """
    Trim the PDB based on the length and the chains specified
    """

    fn = "temp2"
    if not os.path.exists(fn):
        run_command("mkdir "+fn)

    p = PDBParser()
    s = p.get_structure('X', pdbfile)

    outfile = fn+'/'+pdbid+".pdb"
    io = PDBIO()

    class ChainSelector(Select):

        def __init__(self, chains):
            self.chains = chains

        def accept_chain(self, chain):
            return (chain.get_id() in self.chains)

    chains_to_retain = []
    chains_to_retain.append(mhc_chain)
    chains_to_retain.append(pep_chain)
    io.set_structure(s)
    io.save(outfile, ChainSelector(chains_to_retain))

    mhc_chain_obj = s[0]
    pep_chain_obj = s[0]
    for model in s:
        for c in model:
            if mhc_chain == c.id:
                mhc_chain_obj = c
            elif pep_chain == c.id:
                pep_chain_obj = c

    resid_for_chains = get_start_resid_for_chains(pdbfile)

    new_mhc_length = resid_for_chains[mhc_chain] + MHC_TRIM_LENGTH - 1 - len(residues_to_insert)
    new_pep_length = resid_for_chains[pep_chain] + pep_length

    s = p.get_structure('X', outfile)

    class ResSelect(Select):
        def accept_residue(self, res):

            if res.id[1] in mhc_chain_obj and res.id[1] < new_mhc_length:
                return True
            elif res.id[1] in pep_chain_obj and res.id[1] < new_pep_length:
                return True
            else:
                return False

    io.set_structure(s)
    io.save(outfile, ResSelect())

    reordered_file = reorder_chains(pdbid, outfile, mhc_chain, pep_chain)

    return reordered_file


def reorder_chains(filename, infile, mhc_chain, pep_chain):
    """
    Reorder PDB chains such that MHC chain appears first and the peptide chain appears second
    """

    coordinates = defaultdict(list)
    outfile = filename+'_reordered.pdb'

    readfilehandler = open(infile, 'r')
    for line in readfilehandler:
        fields = line.split()
        if len(fields) >= 5:
            chain = fields[21]
            coordinates[chain].append(line)
    readfilehandler.close()

    writefilehandler = open(outfile, 'w')
    for line in coordinates[mhc_chain]:
        writefilehandler.write(line)
    for line in coordinates[pep_chain]:
        writefilehandler.write(line)
    writefilehandler.write("END")
    writefilehandler.close()

    return outfile


def remove_H(pdbfile):

    cmd.load(pdbfile)
    cmd.remove("hydro")
    cmd.save(pdbfile)
    cmd.delete("all")


def rename_chains(pdbfile, mhc_chain, pep_chain):
    """
    Rename PDB chains such that MHC chain is chain A and the peptide chain is chain B (see constants.py)
    """

    cmd.load(pdbfile)
    if mhc_chain != f'{pep_chain}123':
        cmd.alter(f'chain {pep_chain}', f"chain = '{pep_chain}123'")
        cmd.alter(f'chain {mhc_chain}', f"chain = '{DEFAULT_MHC_CHAIN}'")
        cmd.alter(f'chain {pep_chain}123', f"chain = '{DEFAULT_PEP_CHAIN}'")
        cmd.save(pdbfile)
        cmd.delete("all")

    return 'A', 'B'


def update_PDBS_fasta(pdbid, allele, mhc_chain, pep_chain, mhc_seq, pep_seq, outfile, pdbid_date, resolution):
    """
    Update the fasta file to reflect the sequences of the newly downloaded PDB structures
    """

    a = allele.replace("*", "").replace(":", "")
    mhc_header = ">"+pdbid+"|"+mhc_chain+"|"+a+"|"+pdbid_date[pdbid]+"|"+resolution
    pep_header = ">"+pdbid+"|"+pep_chain+"|"+a+"|"+pdbid_date[pdbid]+"|"+resolution

    output = mhc_header+"\n"+mhc_seq+"\n"+pep_header+"\n"+pep_seq+"\n"

    if os.path.isfile(outfile):
        outputfilehandler = open(outfile, 'r')
        if output in outputfilehandler.read():
            print(f"{pdbid} exists in {outfile}")
            outputfilehandler.close()
        else:
            outputfilehandler.close()
            outputfilehandler = open(outfile, 'a')
            outputfilehandler.write(output)
            outputfilehandler.close()

    else:
        outputfilehandler = open(outfile, 'a')
        outputfilehandler.write(output)
        outputfilehandler.close()


def get_statistics(stats, f):
    """
    Get allele and peptide-wise statistics for reporting the new structures downloaded
    """

    allele_dict = {}
    pep_length_dict = {}
    for s in stats:
        pdbid = s[0]
        allele = s[1]
        pep_length = s[2]

        print(pdbid, allele, pep_length, file=f)

        if allele in allele_dict:
            allele_dict[allele] += 1
        else:
            allele_dict[allele] = 1

        if pep_length in pep_length_dict:
            pep_length_dict[pep_length] += 1
        else:
            pep_length_dict[pep_length] = 1

    print(allele_dict, file=f)
    print(pep_length_dict, file=f)


def print_to_file(missing_density_mhc_list, zero_occupancy_mhc_list):
    '''
    Prints to missing density occupancy file for manual curation
    '''

    outputfilehandler = open("PDBS_missing_density_or_zero_occ.txt", "w")
    for mhc in missing_density_mhc_list:
        outputfilehandler.write(mhc+" Missing density in MHC heavy chain\n")

    for mhc in zero_occupancy_mhc_list:
        outputfilehandler.write(mhc+" Zero occupancy found in MHC heavy chain\n")
    outputfilehandler.close()


def cleanup_and_move(curated_files_to_be_removed):
    """
    Remove unwanted structures and directories and put the final PDB files in TRIMMED_TEMPLATES
    and update PDBS.fasta with sequence information
    """

    for pdb in curated_files_to_be_removed:
        run_command("rm -v "+pdb)

    library_dir = 'lib'
    final_templates = library_dir+'/TRIMMED_TEMPLATES'

    if not os.path.exists(library_dir):
        run_command("mkdir "+library_dir)

    if not os.path.exists(final_templates):
        run_command("mkdir "+final_templates)

    run_command("rm rcsb_search_temp.json")
    run_command("rm -r temp1")
    run_command("rm -r temp2")
    run_command("mv *_reordered.pdb lib/TRIMMED_TEMPLATES/")
    run_command("rm *.pdb")
    run_command("rm -r obsolete")
    run_command("rm test*.remodel")
    run_command("rm score.sc")
    run_command("cat PDBS_update.fasta >> lib/PDBS.fasta")
    run_command("rm PDBS_update.fasta")
    run_command("mv PDBS_missing_density_or_zero_occ.txt lib/")


def update_db(rosettainstalldir, start_fasta, date=None, debug_pdbid=None):
    """
    Main method that calls all the other methods
    """

    # gets HLA sequences for HLA typing from HLA_sequences_180.py
    seq2allele = get_hla_sequences()

    # Gets the PDB IDs and respective year
    pdbid_date = fetch_dates_of_pdbs(f"{DATABASE_SCRIPT_PATH}/create_db/rcsb_search.json", debug_pdbid)

    # If no date is given, it defaults to current year
    update_json(f"{DATABASE_SCRIPT_PATH}/create_db/rcsb_search.json", date)

    # Gets PDB IDs as a list
    pdbids = fetch_pdbids(f"{DATABASE_SCRIPT_PATH}/create_db/rcsb_search.json", date, debug_pdbid)

    if pdbids == []:
        print(f"There are no structures >= January 1st, {date}")
        quit()

    # CAUTION: This script is not able to update from an existing fasta file. 
    # All updating is done through a 'restart' of the database so an empty fasta file is provided.
    (pdbs_to_curate, pep_seq_allele_pair) = get_unique_pdbids(start_fasta, pdbids)

    stats_dict = []
    duplicates = defaultdict(list)
    pdbs_curated_and_accepted = {}
    curated_files_to_be_removed = []

    mhc_seq_truncated_list = []
    pdb_has_receptor_list = []
    pdb_not_hla_list = []
    missing_mhc_pep_chain_list = []
    missing_density_pep_list = []
    missing_density_mhc_list = []
    allele_type_not_found_in_seq_list = []
    poor_resolution_list = []
    errored_pdbs_list = []
    pep_over_10AA_list = []
    ignore_receptors_list = []
    zero_occupancy_pep_list = []
    zero_occupancy_mhc_list = []
    chain_truncated = []
    duplicates_num = 0
    approved = 0

    examined_pdbids = []
    print(f'PDBs to curate: {len(pdbs_to_curate)}')

    for pdbid in pdbs_to_curate:

        examined_pdbids.append(pdbid)
        pdbfile, resolution = get_pdb(pdbid)

        try:
            p = PDBParser()
            s = p.get_structure('X', pdbfile)
        except:
            print("PDB reading errored out")
            errored_pdbs_list.append(pdbid)
            continue

        (pdb_has_receptor, pdb_not_hla, mhc_chain, pep_chain, pep_length, remodel, residues_to_insert, pep_over_10AA, b2m_chains, ignore_receptors, mhc_seq_truncated) = align_chains_to_reference(pdbfile)

        if pdb_not_hla:
            pdb_not_hla_list.append(pdbid)

        elif mhc_seq_truncated:
            mhc_seq_truncated_list.append(pdbid)

        elif (pdb_has_receptor and ignore_receptors):
            ignore_receptors_list.append(pdbid)

        elif (pdb_has_receptor and not ignore_receptors):
            pdb_has_receptor_list.append(pdbid)

        elif pep_over_10AA:
            pep_over_10AA_list.append(pdbid)

        elif len(mhc_chain) == 0 or len(pep_chain) == 0:
            missing_mhc_pep_chain_list.append(pdbid)

        if not mhc_seq_truncated and not pdb_not_hla and ((pdb_has_receptor and ignore_receptors) or not pdb_has_receptor) and not pep_over_10AA and len(mhc_chain) > 0 and len(pep_chain) > 0:

            reordered_file = trim_and_save_mhc_pep_chains(pdbfile, pdbid, mhc_chain, pep_chain, pep_length, residues_to_insert)

            if is_missing_bb_heavy_atoms(reordered_file, pep_chain):
                missing_density_pep_list.append(pdbid)
                run_command("rm "+reordered_file)
            elif is_missing_bb_heavy_atoms(reordered_file, mhc_chain):
                missing_density_mhc_list.append(pdbid)
                run_command("rm "+reordered_file)
            elif has_zero_occupancy_atoms(reordered_file, pep_chain):
                zero_occupancy_pep_list.append(pdbid)
                run_command("rm "+reordered_file)
            elif has_zero_occupancy_atoms(reordered_file, mhc_chain):
                zero_occupancy_mhc_list.append(pdbid)
                run_command("rm "+reordered_file)

            else:

                if resolution_of_bb(reordered_file, pep_chain, pep_length):

                    (first_allele, raw_allele) = is_seq_present(seq2allele, get_sequence(reordered_file, mhc_chain))

                    if first_allele != '':

                        if remodel:
                            run_remodel(reordered_file, mhc_chain, pep_chain, residues_to_insert, rosettainstalldir)
                        else:
                            renumber_residues(reordered_file)

                        mhc_chain, pep_chain = rename_chains(reordered_file, mhc_chain, pep_chain)

                        allele_seq = get_sequence(reordered_file, mhc_chain)
                        pep_seq = get_sequence(reordered_file, pep_chain)
                        pMHC = pep_seq+"_"+first_allele

                        if (pMHC not in pdbs_curated_and_accepted) and (pMHC not in pep_seq_allele_pair):
                            approved += 1
                            pdbs_curated_and_accepted[pMHC] = pdbid
                            update_PDBS_fasta(pdbid, raw_allele, mhc_chain, pep_chain, allele_seq, pep_seq, "PDBS_update.fasta", pdbid_date, resolution)
                            stats_dict.append((pdbid, raw_allele, pep_length))

                        else:
                            curated_files_to_be_removed.append(reordered_file)
                            if pdbid not in duplicates[pMHC]:
                                duplicates[pMHC].append(pdbid)
                                duplicates_num += 1
                            if pdbs_curated_and_accepted[pMHC] not in duplicates[pMHC]:
                                duplicates[pMHC].append(pdbs_curated_and_accepted[pMHC])

                    else:
                        curated_files_to_be_removed.append(reordered_file)
                        allele_type_not_found_in_seq_list.append(pdbid)
                else:
                    curated_files_to_be_removed.append(reordered_file)
                    poor_resolution_list.append(pdbid)


    print_to_file(missing_density_mhc_list, zero_occupancy_mhc_list)
    cleanup_and_move(curated_files_to_be_removed)

    with open('lib/output_stats.txt', 'w') as f:
        get_statistics(stats_dict, f)
        print("Total PDBs started with:", len(pdbs_to_curate), file=f)
        print("PDBs is not a HLA:", len(pdb_not_hla_list), file=f)
        print(','.join(pdb_not_hla_list), file=f)
        print("PDBS have short heavy chain:", len(mhc_seq_truncated_list), file=f)
        print(",".join(mhc_seq_truncated_list), file=f)
        print("PDBs with peptides over 10 AA:", len(pep_over_10AA_list), file=f)
        print(','.join(pep_over_10AA_list), file=f)
        print("PDBs have receptor closer to peptide:", len(pdb_has_receptor_list), file=f)
        print(','.join(pdb_has_receptor_list), file=f)
        print("PDBs with receptors away from the peptide (double counted):", len(ignore_receptors_list), file=f)
        print(','.join(ignore_receptors_list), file=f)
        print("PDBs have missing MHC or peptide chain:", len(missing_mhc_pep_chain_list), file=f)
        print(','.join(missing_mhc_pep_chain_list), file=f)
        print("PDBs have missing density in peptide:", len(missing_density_pep_list), file=f)
        print(','.join(missing_density_pep_list), file=f)
        print("PDBs have missing density in heavy chain:", len(missing_density_mhc_list), file=f)
        print(','.join(missing_density_mhc_list), file=f)
        print("PDBs have zero occupancy atoms in peptide:", len(zero_occupancy_pep_list), file=f)
        print(','.join(zero_occupancy_pep_list), file=f)
        print("PDBs have zero occupancy atoms in heavy chain:", len(zero_occupancy_mhc_list), file=f)
        print(','.join(zero_occupancy_mhc_list), file=f)
        print("PDBs allele type not found in the sequence database:", len(allele_type_not_found_in_seq_list), file=f)
        print(','.join(allele_type_not_found_in_seq_list), file=f)
        print("Poor resolution PDBs (double bb):", len(poor_resolution_list), file=f)
        print(','.join(poor_resolution_list), file=f)
        print("Duplicate PDBs:", file=f)
        dups = 0
        for key in duplicates:
            print(key, '->', ','.join(duplicates[key]), file=f)
            dups += len(duplicates[key]) - 1
        print("Total duplicates:", dups, file=f)
        print("PDBs errored due to some issue while downloading:", len(errored_pdbs_list), file=f)
        print(','.join(errored_pdbs_list), file=f)

        executionTime = (time.time() - startTime)
        print('Script ran in', str(round(executionTime, 3)), 'seconds', file=f)

if __name__ == "__main__":

    args = get_args()
    update_db(args.rosettainstalldir, args.start_fasta, args.date, args.debug_pdbid)
