#!/bin/bash


mass_defs=("/mean" "/crit")
cost_types=("" "/abund")
fit_types=("/rd3.param" "/rd_rm3.param" "/rm3.param")
complete=()
incomplete=()
for mass_def in "${mass_defs[@]}"; 
do
    for cost_type in "${cost_types[@]}"; do
    	for fit_type in "${fit_types[@]}"; do
	    param_file="${1}${mass_def}${cost_type}${fit_type}"
    	    ./core_fit "${param_file}"
	    exit_code=$?
	    echo "\n\nexitcode: $?"
	    if [ "$?" -eq "0" ] 
	    then
		complete+=($param_file)
		./calc_likelihood_bounds.py "${param_file}"
	    else
		incomplete+=($param_file)
	    fi
	    echo "complete: ${complete[@]}"
	    echo "incomplete: ${incomplete[@]}"
    	done
    done
done
