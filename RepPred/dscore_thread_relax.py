import sys
import pandas as pd
import numpy as np
import pyrosetta
from pyrosetta import pose_from_pdb
from pyrosetta import *


def dscore_calc(angle1, angle2):
	distance = 2.0 * np.sum(1.0 - np.cos(np.radians(angle1) - np.radians(angle2)))
	return distance

def calc_dihedrals(pdbfile):

	cinfo = ("B", 4, 7) # peptide is chain B and middle of peptide is P4 to P7

	angle1 = []	
	pose = pose_from_pdb(pdbfile)
	p1 = pose.pdb_info().pdb2pose(cinfo[0], cinfo[1])
	p0 = pose.pdb_info().pdb2pose(cinfo[0], cinfo[2])
	for resi_num in range(p1,p0+1):
		angle1 += [pose.phi(resi_num), pose.psi(resi_num)]
	return np.array(angle1)
import sys

pyrosetta.init()
data=pd.DataFrame(columns=['dscore','target_on_template'])

targname=sys.argv[1]   
tempname=sys.argv[2]   
num=sys.argv[3]   

targfile='./run-dir/relaxed_structures/'+targname+'_on_'+tempname+'/'+targname+'_on_'+tempname+'_reordered_relaxed_000'+str(num)+'.pdb'
tempfile='./run-dir/'+targname+'/'+targname+'_on_'+tempname+'_reordered.pdb'
tar_angle = calc_dihedrals(targfile)
temp_angle = calc_dihedrals(tempfile)
dscore= dscore_calc(tar_angle, temp_angle)
print(targname,tempname)
data=data.append({'target_on_template':f'{targname}_on_{tempname}_{num}','dscore':dscore}, ignore_index=True)
data.to_csv(f"dscores_thread_relax/dscore_{targname}_{tempname}_{num}.csv",index=False)




