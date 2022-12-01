'''
score_decoys.py
PyRosetta4/python3 script for generating score function features for SVR selection.
Works for any length of peptide, but note that the selection models only predict for nonameric peptides.
-- GLJ Keller
'''
import os
import pyrosetta
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.scoring import *
from constants import *

pyrosetta.init(extra_options = "-extrachi_cutoff 12 -ex1 -ex2 -ex3 -ignore_zero_occupancy false")
scorefxn = pyrosetta.create_score_function('ref2015_cart')

def scoring(model, name, outputfile):

	outData = feature_headers(9)

	decoy = pose_from_pdb(model)

	p1 = decoy.pdb_info().pdb2pose(PEP_CHAIN, 1)
	pO = decoy.pdb_info().pdb2pose(PEP_CHAIN, 9)

	score_data = []
	score = scorefxn(decoy)
	energies = decoy.energies()
	residue_sasa = pyrosetta.Vector1([1.0, 2.0])
	residue_hydro_sasa = pyrosetta.Vector1([3.0, 4.0])
	score_data.extend([str(decoy.energies().residue_total_energies(a)[fa_atr]) for a in range(p1, pO+1)])
	score_data.extend([str(decoy.energies().residue_total_energies(b)[fa_rep]) for b in range(p1, pO+1)])
	score_data.extend([str(decoy.energies().residue_total_energies(g)[fa_sol]) for g in range(p1, pO+1)])
	score_data.extend([str(decoy.energies().residue_total_energies(d)[fa_intra_sol_xover4]) for d in range(p1, pO+1)])
	score_data.extend([str(decoy.energies().residue_total_energies(e)[lk_ball_wtd]) for e in range(p1, pO+1)])
	score_data.extend([str(decoy.energies().residue_total_energies(j)[fa_intra_rep]) for j in range(p1, pO+1)])
	score_data.extend([str(decoy.energies().residue_total_energies(e_again)[fa_elec]) for e_again in range(p1, pO+1)])
	score_data.extend([str(decoy.energies().residue_total_energies(k)[hbond_bb_sc]) for k in range(p1, pO+1)])
	score_data.extend([str(decoy.energies().residue_total_energies(i)[hbond_sc]) for i in range(p1, pO)])
	score_data.extend([str(decoy.energies().residue_total_energies(y)[rama_prepro]) for y in range(p1+1, pO)])
	score_data.extend([str(decoy.energies().residue_total_energies(yy)[omega]) for yy in range(p1, pO)])
	score_data.extend([str(decoy.energies().residue_total_energies(yyz)[p_aa_pp]) for yyz in range(p1+1, pO)])
	score_data.extend([str(decoy.energies().residue_total_energies(o)[fa_dun]) for o in range(p1, pO+1)])
	pyrosetta.rosetta.core.scoring.calc_per_res_hydrophobic_sasa(decoy, residue_sasa, residue_hydro_sasa, 1.4, False)
	score_data.extend([str(residue_sasa[s]) for s in range(p1, pO+1)])
	score_data.extend([str(residue_hydro_sasa[h]) for h in range(p1, pO+1)])

	data = [name, str(score)]
	data.extend(score_data)
	outData.append(','.join(data))

	if os.path.exists(outputfile):
		with open(outputfile, 'a') as OF:
			print(''.join(outData[1]), file=OF)
	else:
		with open(outputfile, 'a') as OF:
			print('\n'.join(outData), file=OF)

def feature_headers(end_res):
	outData = [['target_on_template','totalScore']]
	outData[0].extend([f'fa_atr{a}' for a in range(1, end_res+1)]) # Lennard-Jones attractive between atoms in different residues
	outData[0].extend([f'fa_rep{b}' for b in range(1, end_res+1)]) # Lennard-Jones repulsive between atoms in different residues
	outData[0].extend([f'fa_sol{g}' for g in range(1, end_res+1)]) # Lazaridis-Karplus solvation energy
	outData[0].extend([f'fa_intra_sol_xover4{d}' for d in range(1, end_res+1)]) # Intra-residue Lazaridis-Karplus solvation energy
	outData[0].extend([f'lk_ball_wtd{e}' for e in range(1, end_res+1)]) # Asymmetric solvation energy
	outData[0].extend([f'fa_intra_rep{j}' for j in range(1, end_res+1)]) # Lennard-Jones repulsive between atoms in the same residue
	outData[0].extend([f'fa_elec{e_again}' for e_again in range(1, end_res+1)]) # Coulombic electrostatic potential with a distance-dependent dielectric
	outData[0].extend([f'hbond_bb_sc{k}' for k in range(1, end_res+1)]) # Sidechain-backbone hydrogen bond energy
	outData[0].extend([f'hbond_sc{i}' for i in range(1, end_res)]) # Sidechain-sidechain hydrogen bond energy
	outData[0].extend([f'rama_prepro{y}' for y in range(2, end_res)]) # Ramachandran preferences (with separate lookup tables for pre-proline positions and other positions)
	outData[0].extend([f'omega{yy}' for yy in range(1, end_res)]) # Omega dihedral in the backbone. A Harmonic constraint on planarity with standard deviation of ~6 deg.
	outData[0].extend([f'p_aa_pp{yyz}' for yyz in range(2, end_res)]) # Probability of amino acid, given torsion values for phi and psi
	outData[0].extend([f'fa_dun{o}' for o in range(1, end_res+1)]) # Internal energy of sidechain rotamers as derived from Dunbrack's statistics
	outData[0].extend([f'sasa{s}' for s in range(1, end_res+1)]) # Solvent-accessible surface area of residue based on 1.4A probe.
	outData[0].extend([f'hydro_sasa{h}' for h in range(1, end_res+1)]) # Solvent-accessible surface area of hydrophobic atoms in residue based on 1.4A probe.
	outData = [','.join(outData[0])]
	return outData
