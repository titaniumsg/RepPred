#       Sgourakis Lab
#   Author: Sagar Gupta
#   Date: March 18, 2023
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
9mer peptides with anchor residues P2 and P9 bound to an HLA

'''

def get_args():
   
    '''
    Parse command line arguments
    '''

    parser = argparse.ArgumentParser(description="Structural modeling of peptide/HLA complexes")
    parser.add_argument("-seq_file", help="File containing HLA sequence (180 residues) followed by peptide sequence (9 residues) on a new line", required=True)
    parser.add_argument("-output_dir", help="Output directory (Highly recommend for this to be empty!)", default="example")
    parser.add_argument("-prefix", help="Prefix used for all files", default="test")

    parser.add_argument("-num_iters", help="Number of FastRelax rounds to run (at least 3 is recommended)", default="4", type=str)

    args = parser.parse_args()

    return args

def main(seq_file, output_dir, prefix, num_iterations):
    
    def read_seq_file(seq_file):
        with open(seq_file, "r") as openfile:
            sequences = openfile.readlines()
            return sequences[0].strip(), sequences[1].strip()

    hla_sequence, peptide_sequence = read_seq_file(seq_file)

    assert len(peptide_sequence) == 9, f"Peptide must be length 9. A peptide of length {len(peptide_sequence)} is not supported"
    assert len(hla_sequence) == 180, f"HLA must be length 180. An HLA of length {len(hla_sequence)} is not supported"
    
    print(f"Predicting the structure of the following peptide/HLA sequence:\n{hla_sequence}\n{peptide_sequence}\n")

    model = pickle.load(open(f'{REP_PRED_DIR}/regression/model.pkl', 'rb'))
    scaler = pickle.load(open(f'{REP_PRED_DIR}/regression/scaler.pkl', 'rb'))
    predicted_dscore_dict = {}

    structural_templates = sorted(glob.glob(f"{REP_PRED_DIR}/templates/*.pdb"))
    for template in tqdm(structural_templates):
        
        threaded_filename = thread(template, f"{output_dir}/{prefix}", peptide_sequence, hla_sequence)

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
    
    run_command(f"mkdir -p {output_dir}/workdir")

    run_command(f"rm {output_dir}/test.grishin")

    run_command(f"mv {output_dir}/*.log {output_dir}/workdir")
    run_command(f"mv {output_dir}/*.pdb {output_dir}/workdir")
    run_command(f"mv {output_dir}/*.csv {output_dir}/workdir")
    run_command(f"mv {output_dir}/*.sc {output_dir}/workdir")

    best_model_base = os.path.basename(best_model)

    run_command(f"cp {output_dir}/workdir/{best_model_base} {output_dir}")
    
    print(f"The best model of {output_dir} has a predicted D-score of {predicted_dscore_dict[best_model]:.2f}\nThe modeled PDB file can be found at {output_dir}/{best_model_base}")

if __name__ == "__main__":

    args = get_args()
    
    main(args.seq_file, args.output_dir, args.prefix, args.num_iters)
