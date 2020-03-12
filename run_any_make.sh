#!/bin/bash
echo $#
if [ $# -l 1 ]; then
    echo 'run mode: ./run_any_make.sh exe make_param_sh'
    exit
fi
exe=$1
make_script=$2
if [ -z "$3"]; then
    zoom="false"
else
    zoom=$3
fi
results=`bash ${make_script} false`


for result in ${results}; do
    if [ $zoom == "zoom" ]; then
	result="${result/.param/_zoom.param}"
    fi
    echo ""
    echo ""
    echo ""
    echo ${result}
    nice -n 13 ${exe} ${result}
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
