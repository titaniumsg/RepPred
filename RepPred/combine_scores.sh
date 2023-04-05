
mkdir -p model_predictions

head -1 score_header.csv>>scores/model_scores.csv
for targets in input_sequence/*.txt; do
	targname=$(echo ${targets} | cut -f 2 -d '/' | cut -f 1 -d '_')
	for templates in template_pdbs/*.pdb; do
		tempname=$(echo ${templates} | cut -f 2 -d '/' | cut -f 1 -d '_')
		echo ${targname} ${tempname}
		tail -1 scores/score_${targname}_on_${tempname}.csv >>scores/model_scores.csv
	done
done
#python svr_predictions.py
