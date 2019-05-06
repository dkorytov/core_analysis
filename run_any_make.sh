#!/bin/bash

exe=$1
make_script=$2

results=`bash ${make_script} false`

for result in ${results}; do
    echo ${result}
    ${exe} ${result}
done
# # param_pattern=`sed
# for mass in "mean" "crit";do
#     for model in "rd" "rm" "rd_rm"; do
# 	for fit in "" ; do #"/abund"
# 	    param=`echo ${param_pattern} | sed "s/@mass@/${mass}/g" | sed "s/@model@/${model}/g" | sed "s%@fit@%${fit}%g"`
# 	    echo ${param}
# 	    $exe ${param}
# 	done
#     done
# done
