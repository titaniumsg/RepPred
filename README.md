# HLA3DB and RepPred

 #### For a detailed description of the methods, please refer to the article:
    "Gupta, et al. HLA3DB: comprehensive annotation of peptide/HLA complexes enables blind structure prediction of T cell epitopes. 
    bioRxiv. doi.org/10.1101/2023.03.20.533510"
Corresponding author: Nikolaos G. Sgourakis, nikolaos.sgourakis@pennmedicine.upenn.edu


#### System requirements:
    - MacOS/Linux (tested on MacOS 12.6.3 and RHEL 8.5)
    - Python 3.8+ (tested on 3.8.15)
    - Rosetta (tested on 2020.08) – compiled with MPI is recommended
    - Anaconda (tested on 4.11.0)
        - pip
            - biopython
            - scikit-learn
        - pandas
        - platformdirs
        - pymol
        - tqdm
        - pyrosetta 

HLA3DB can be downloaded from https://hla3db.research.chop.edu/ or generated via a normal desktop computer.

To run RepPred, we recommend using an HPC. We have configured our code for a slurm scheduler.


#### Installation:
1.  Modify the `environment.yml` file to include the appropriate `username` and `password` for PyRosetta installation. 
A license can be obtained via: https://els2.comotion.uw.edu/product/pyrosetta
2. Create a conda environment using the provided YAML file using the following command:
`conda env create -f environment.yml`
This step takes roughly 5 minutes on a normal desktop computer, depending on internet speed.
3. Obtain Rosetta via https://www.rosettacommons.org/software/license-and-download



### HLA3DB

An auto-updating, curated version of HLA3DB is readily accessible via https://hla3db.research.chop.edu/ for download.
We also provide the option to locally generate the database below.

The database provided in HLA3DB/create_db/lib/ was updated as of April 4th, 2023.
The corresponding statistics can be found in db_stats/database.csv

Please ensure the appropriate conda environment is activated:
    `conda activate hla3db`

#### DATABASE CREATION:

1. Modify `create_db.sh` by changing `/your_rosetta_dir/` to your local Rosetta folder
2. Run `bash create_db.sh`. Script takes ~60 mins to run on a normal desktop computer.
3. Curated structures can be found in `create_db/lib/TRIMMED_TEMPLATES/`.

#### HLA3DB STATISTICS:
Must be done *AFTER* creating HLA3DB
1. Run `bash db_stats.sh`. Script takes ~30 mins to run on a normal desktop computer.
2. Output can be found in `db_stats/database.csv`.



### RepPred

Please ensure the appropriate conda environment is activated:
    `conda activate hla3db`

An example input can be found in the `input_sequence` folder and
an example output can be found in the `model_predictions` folder.

#### INPUT

The required input from the user is the sequence of HLA-I (180 residues) and peptide (9 residues) in a .txt format 
(Format provided in `input_sequence` folder)

#### CODE USAGE

The code can be run in two steps:

Step 1:

   `bash run_RepPred.sh <option> <rosetta_dir>`

   Code is written to run either on a hpc-cluster with slurm scheduler (`<option>`: hpc) 
   or a local machine (`<option>`: local)
The local Rosetta directory also need to be provided via `<rosetta_dir>`

Step 2:    

   Proceed to this step *ONLY* after step 1 is completed.

`python svr_predictions.py`

   The models will be generated in the `model_predictions` folder.
   This folder also contains information about the predicted dscore and the template 
   used to produce the model in .csv format

#### RUNTIME
	
We *strongly* reccomend running the code in an hpc environment, this should take approximately 
   20-30 minutes depending on the query peptide. 

When run on a local machine the runtime will be approximately 15 hours.

#### NOTE
     
In the modeling protocol, if the threaded model generated from a template is structurally
dissimilar from  the template (D-score > 1.5) that model is not included for further process.
Hence a message like "cannot open 'scores/score_7RE8_on_5ENW.csv' for reading " is normal to appear.
