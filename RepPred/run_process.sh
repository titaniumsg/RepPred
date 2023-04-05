#!/bin/bash
start_time=`date +%s`
target=$1
template=$2
rundir=$3
relax_method=$4
cluster_opt=$5
rosetta_dir=$6
python run_thread.py -run_dir ${rundir} -structure_dir template_pdbs -target_pdbid ${target} -template_pdbid ${template}

pdbname=${rundir}/${target}/${target}_on_${template}_reordered.pdb
outdir=$(echo ${rundir}/relaxed_structures/${target}_on_${template})

if [[ "${cluster_opt}" == "hpc" ]]; then
	mpirun -np 4 ${rosetta_dir}/source/bin/relax.mpi.linuxgccrelease -s ${pdbname} -database ${rosetta_dir}/database -score:weights ref2015_cart  -relax:cartesian  -relax:constrain_relax_to_start_coords -relax:ramp_constraints false -out:suffix _relaxed -nstruct 4 -out:path:score ${outdir} -out:path:pdb ${outdir}

else
	${rosetta_dir}/source/bin/relax.macosclangrelease -s ${pdbname} -database ${rosetta_dir}/database -score:weights ref2015_cart  -relax:cartesian  -relax:constrain_relax_to_start_coords -relax:ramp_constraints false -out:suffix _relaxed -nstruct 4 -out:path:score ${outdir} -out:path:pdb ${outdir}

fi

keep_struct=0     # best structure is the one with lowest energy and less than cutoff dscore from template
for i in {1..4};do                 # check for all relaxes
	echo ${target} ${template} ${i}
	if (( keep_struct == 0 ));then   #***IMPORTANT****change the line below if it is not cartesian relax
		best=$(grep ${target}_on_${template} ${outdir}/score_relaxed.sc | awk '{print $23" "$2-$4}' | sort -n -k 2 | head -n ${i}| tail -n 1 | awk '{print $1}')    # get best energy structure
		echo 'BEST STRUCTURE',${best}
		
		relaxnum=$(($(echo ${best} | cut -f 6 -d '_'))) # get corresponding relax number
		bestfile=$(echo ${outdir}/$best.pdb)
		python dscore_thread_relax.py ${target} ${template} ${relaxnum}
		dscore_val=$(tail -1 dscores_thread_relax/dscore_${target}_${template}_${relaxnum}.csv| cut -f 1 -d ',')
		echo 'DSCORE',${dscore_val}
		cutoff=1.5
		if (( $(echo "${dscore_val} < ${cutoff}" |bc -l ) )); then
 			keep_struct=1
 			echo 'KEEP THIS STRUCTURE'
 			echo ${target}_${template},${dscore_val},$i >>dscores_thread_relax/selected_dscores.csv
 			newfile=$(echo ${rundir}/relaxed_structures/best_score/${target}_on_${template}.pdb)
 			cp -vf $bestfile ${newfile}
		fi
 	fi
done

if (( keep_struct == 1 ));then
	python per_residue_score_single.py $target $template
fi
end_time=`date +%s`

echo Execution_time_in seconds, `expr $end_time - $start_time`>>time_data/time_${target}_${template}.csv 
