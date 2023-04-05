import sys,os
import subprocess

seq_path=sys.argv[1]
sys_status=sys.argv[2]
print("seq path :",seq_path)

if not  os.path.exists(seq_path):
    print("Given sequence path does not exist.")
    sys.exit()

#make necessary folders
dirs=subprocess.run(f"bash run_rosetta_model.sh {seq_path} {sys_status}",  shell=True)


