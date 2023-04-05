#       Sgourakis Lab
#   Author: Sagar Gupta
#   Date: April 28, 2022
#   Email: sagarg@sas.upenn.edu

# import custom libraries
from misc.constants import *
from create_db.update_db import update_db
from create_db.manual_addition.add_to_db import add_to_db

import time
startTime = time.time()

def main():

    update_db(rosettainstalldir=ROSETTA_INSTALL_DIR, start_fasta=f'{DATABASE_SCRIPT_PATH}/create_db/PDBS_empty.fasta')

    add_to_db()

    executionTime = (time.time() - startTime)
    print('Script ran in', str(round(executionTime, 3)), 'seconds')

if __name__ == "__main__":

    main()
