#!/bin/bash

make_runscript(){
rundir=$1
target=$2
template=$3
relax_method=$4
cluster_opt=$5
rosetta_dir=$6
   
cat <<EOF > jobfiles/submit_script-${target}-${template}.sh
#!/bin/sh
#
#SBATCH --job-name=job-${target}-${template}
#SBATCH --output=jobfiles/thread-${target}-${template}.txt
#SBATCH --error=jobfiles/error-${target}-${template}.txt
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=4
module load openmpi/gcc/
bash run_process.sh ${target} ${template} ${rundir} ${relax_method} ${cluster_opt} ${rosetta_dir}
exit
EOF
}

seqpath=${1}
cluster_opt=${2}
rosetta_dir=${3}

mkdir -p run-dir
mkdir -p run-dir/relaxed_structures
mkdir -p run-dir/relaxed_structures/best_score
mkdir -p scores
mkdir -p model_predictions
mkdir -p jobfiles
mkdir -p dscores_thread_relax
mkdir -p time_data


relax_method='cartesian'

for targets in ${seqpath}/*.txt; do
	targname=$(echo ${targets} | cut -f 2 -d '/' | cut -f 1 -d '_')
	mkdir run-dir/${targname}
	for templates in template_pdbs/*.pdb; do
		tempname=$(echo ${templates} | cut -f 2 -d '/' | cut -f 1 -d '_')
		mkdir run-dir/relaxed_structures/${targname}_on_${tempname}
		if [ ${cluster_opt} == 'hpc' ]
		then
			echo "Selected option is hpc."
        	make_runscript run-dir ${targname} ${tempname} ${relax_method} ${cluster_opt} ${rosetta_dir}
	    	sbatch jobfiles/submit_script-${target}-${template}.sh
		else
			echo "Selected option is local machine... this may take longer than hpc option"
			bash run_process.sh ${targname} ${tempname} run-dir ${relax_method} ${cluster_opt} ${rosetta_dir}
		fi
	done
done
