#!/bin/bash

exe=$1
param_pattern=`echo $2 | sed "s%mean%@mass@@fit@%g" | sed "s/rd/@model@/g"`

# param_pattern=`sed
for mass in "mean" "crit";do
    for model in "rd" "rm" "rd_rm"; do
	for fit in "" "/abund"; do 
	    param=`echo ${param_pattern} | sed "s/@mass@/${mass}/g" | sed "s/@model@/${model}/g" | sed "s%@fit@%${fit}%g"`
	    echo ${param}
	    $exe ${param}
	done
    done
done
