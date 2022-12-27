#       Sgourakis Lab
#   Author: Sagar Gupta
#   Date: November 29, 2022
#   Email: sagarg@sas.upenn.edu

# import required libraries
import subprocess
import platform
import os

'''
Relaxes threaded file
'''

def run_command(command, ignore=False):
    '''
    Helps run command line operations
    '''

    cmd = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = cmd.communicate()

    if cmd.returncode != 0:

        if ignore:
            return False

        print("--- FAIL ---\n")
        print(command)
        if output[0] != None:
            print(str(output[0], 'utf-8'))
        if output[1] != None:
            print(str(output[1], 'utf-8'))
        print("\n--- FAIL ---")
        
        return False

    else:
        if output[0] != None:
            output = output[0].decode('utf-8')
        
        return output

def submit_to_cluster(command):

    with open("submit.sh", "w") as bashfile:
        bashfile.write("#!/bin/bash\n")
        bashfile.write("#SBATCH --time=2-23:59:59\n")
        # bashfile.write(". /etc/profile\n")
        # bashfile.write("module load OpenMPI/3.1.4-GCC-8.3.0\n")
        # bashfile.write("unset module\n")

        # bashfile.write(f"mpirun -n 3 {command}")
        bashfile.write(f"{command}")

    run_command("sbatch submit.sh")


def get_extension():
    
    os_system = platform.system()
    extension = ''
    if os_system == 'Linux':
        extension = 'linuxgccrelease'
    elif os_system == 'Darwin':
        extension = 'macosclangrelease'
    
    return extension

def relax(pdb_file, rosetta_dir, num_iterations=4):

    extension = get_extension()
    output_dir = os.path.split(pdb_file)[0]
    pdb_file_base = os.path.basename(pdb_file)
    pdb_file_no_extension = os.path.splitext(pdb_file_base)[0]
    logfile = f"{output_dir}/{pdb_file_no_extension}.log"
    
    print(f"Relaxing {pdb_file} {num_iterations} time(s)\nLog file: {logfile}")
    # run_command(f"{rosetta_dir}/main/source/bin/relax.mpi.{extension} -database {rosetta_dir}/main/database -s {pdb_file} -nstruct {num_iterations} -out:path:score {output_dir} -out:path:pdb {output_dir} -score:weights ref2015_cart -relax:constrain_relax_to_start_coords -relax:cartesian -relax:ramp_constraints false -out:suffix _relaxed > {logfile}")
    submit_to_cluster(f"{rosetta_dir}/main/source/bin/relax.mpi.{extension} -database {rosetta_dir}/main/database -s {pdb_file} -nstruct {num_iterations} -out:path:score {output_dir} -out:path:pdb {output_dir} -score:weights ref2015_cart -relax:constrain_relax_to_start_coords -relax:cartesian -relax:ramp_constraints false -out:suffix _relaxed > {logfile}")
    
    # relaxed_threaded_total_score = {}
    
    # for index in range(1, int(num_iterations)+1):
    #     relaxed_threaded_filename = f"{pdb_file_no_extension}_relaxed_{str(index).zfill(4)}"
    #     total_score = run_command(f"cat {output_dir}/score_relaxed.sc | grep {relaxed_threaded_filename}"+" | tail -1 | awk '{print $2}'")
    #     if total_score == '': continue
    #     total_score = float(total_score.strip())
    #     relaxed_threaded_total_score[relaxed_threaded_filename] = total_score

    # relaxed_threaded_best = min(relaxed_threaded_total_score, key=relaxed_threaded_total_score.get)

    # return f"{output_dir}/{relaxed_threaded_best}.pdb", relaxed_threaded_best
