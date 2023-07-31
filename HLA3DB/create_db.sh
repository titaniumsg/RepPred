#!/bin/bash

# setup
export PYTHONPATH=$PWD
sed -e 's+/database/'+$PWD+'g' -e 's+/rosetta/+/Users/titanium/Downloads/rosetta/+g' misc/constants_template.py > misc/constants2.py
mv misc/constants2.py misc/constants.py

cd create_db
rm -rfv temp1 temp2 obsolete lib
rm -fv *.pdb test*.remodel rcsb_search_temp.json score.sc PDBS_update.fasta
mkdir -p lib/TRIMMED_TEMPLATES

# create the database
python update_db.py
cd manual_addition
python add_to_db.py