Ensure the appropriate conda environment is activated:
    conda activate RepPred

To create HLA3DB: 
1. Modify "create_db.sh" by changing "/your_rosetta_dir/" 
2. Run "bash create_db.sh". Script takes ~60 mins to run on a normal desktop computer.
3. Structures can be found in create_db/lib/TRIMMED_TEMPLATES/

To generate HLA3DB statistics (must be done AFTER creating HLA3DB):
1. Run "bash db_stats.sh". Script takes ~5 mins to run on a normal desktop computer.
2. Output can be found in "db_stats/database.csv"