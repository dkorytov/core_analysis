#!/bin/bash
echo $#
if [ $# -lt 1 ]; then
    echo 'run mode: ./run_any_make.sh exe make_param_sh'
    exit
fi
exe=$1
make_script=$2
if [ -z "$3" ]; then
    zoom="false"
else
    zoom=$3
    if [ -z "$4" ]; then
	zoom2='false'
    else
	zoom2=$4
    fi
fi
results=`bash ${make_script} false`


for result in ${results}; do
    if [ $zoom == "zoom" ]; then
	result="${result/.param/_zoom.param}"
    fi
    if [ $zoom2 == "zoom" ]; then
	result="${result/.param/_zoom.param}"
    fi
    echo ""
    echo ""
    echo ""
    echo ${result}
    nice -n 20 ${exe} ${result}
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
