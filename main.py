#       Sgourakis Lab
#   Author: Sagar Gupta
#   Date: November 29, 2022
#   Email: sagarg@sas.upenn.edu

# import required libraries
import argparse
import glob
import pickle
import pandas as pd
from tqdm import tqdm
import csv
import os

# import custom libraries
from constants import *
from thread import thread
from relax import relax, run_command
from scoring import scoring

'''
main.py is the entry point of RepPred

This method is designed specifically to work with
9mer peptides with anchor residues P2 and P9 bound to HLA-A*02:01

'''

def get_args():
   
    '''
    Parse command line arguments
    '''

    parser = argparse.ArgumentParser(description="Structural modeling of peptide/HLA-A*02:01 complexes")
    parser.add_argument("-peptide_sequence", help="Sequence of desired peptide", required=True)
    parser.add_argument("-output_dir", help="Output directory (Highly recommend for this to be empty!)", default="example")
    parser.add_argument("-prefix", help="Prefix used for all files", default="test")

    parser.add_argument("-num_iterations", help="Number of FastRelax rounds to run (at least 3 is recommended)", default="4", type=str)

    args = parser.parse_args()

    return args

def main(peptide_sequence, output_dir, prefix, num_iterations):
    
    model = pickle.load(open(f'{REP_PRED_DIR}/regression/model.pkl', 'rb'))
    scaler = pickle.load(open(f'{REP_PRED_DIR}/regression/scaler.pkl', 'rb'))
    predicted_dscore_dict = {}

    structural_templates = sorted(glob.glob(f"{REP_PRED_DIR}/templates/*.pdb"))
    for template in tqdm(structural_templates):
        if template != f'{REP_PRED_DIR}/templates/1A1M_reordered_const_cart_relaxed_0001.pdb': continue
        
        threaded_filename = thread(template, f"{output_dir}/{prefix}", peptide_sequence, A0201_HLA_SEQUENCE)

        relaxed_threaded_filename, relaxed_threaded_best = relax(threaded_filename, ROSETTA_DIR, num_iterations)

        scorefile = f"{output_dir}/{relaxed_threaded_best}.csv"

        scoring(relaxed_threaded_filename, relaxed_threaded_best, scorefile)

        data = pd.read_csv(scorefile)
        data.drop(['target_on_template', 'totalScore'], axis=1, inplace=True)
        data_scaled = scaler.transform(data)
        predicted_dscore = model.predict(data_scaled)
        predicted_dscore_dict[relaxed_threaded_filename] = predicted_dscore[0]
    
    predicted_dscore_dict = dict(sorted(predicted_dscore_dict.items(), key=lambda item: item[1])) # sort dictionary by predicted D-score
    best_model = min(predicted_dscore_dict, key=predicted_dscore_dict.get)
    
    with open(f"{output_dir}/{prefix}_pred_dscore.csv", "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['pep_seq', 'template', 'pred_dscore'])
        for template, pred_dscore in predicted_dscore_dict.items():
            pdb_file_base = os.path.basename(template)
            writer.writerow([peptide_sequence, pdb_file_base, pred_dscore])
    
    run_command(f"mkdir {output_dir}/workdir")

    run_command(f"rm {output_dir}/test.grishin")

    run_command(f"mv {output_dir}/*.log {output_dir}/workdir")
    run_command(f"mv {output_dir}/*.pdb {output_dir}/workdir")
    run_command(f"mv {output_dir}/*.csv {output_dir}/workdir")
    run_command(f"mv {output_dir}/*.sc {output_dir}/workdir")

    best_model_base = os.path.basename(best_model)

    run_command(f"cp {output_dir}/workdir/{best_model_base} {output_dir}")
    
    print(f"The best model of {peptide_sequence}/HLA-A*02:01 has a predicted D-score of {predicted_dscore_dict[best_model]:.2f}\nThe modeled PDB file can be found at {output_dir}/{best_model_base}")

if __name__ == "__main__":

    args = get_args()
    assert len(args.peptide_sequence) == 9, f"A peptide of length {len(args.peptide_sequence)} is not supported"
    main(args.peptide_sequence, args.output_dir, args.prefix, args.num_iterations)
