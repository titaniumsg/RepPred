#!/bin/bash

# setup
export PYTHONPATH=$PWD
sed -e 's+/database/'+$PWD+'g' misc/constants_template.py > misc/constants2.py
mv misc/constants2.py misc/constants.py 

cd db_stats
python generate_stats.py