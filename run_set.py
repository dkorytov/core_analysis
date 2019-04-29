#!/bin/bash

# removing display so that it can run in a screen more easily (i.e. I
# don't have to start a new ssh terminal without x-forwarding)
unset DISPLAY

param=$1
param_zoom=`echo ${param} | sed 's/.param/_zoom.param/'`

echo $param
echo $param_zoom

# ./run_any.sh ./core_fit $param
./run_any.sh ./plot_zmrs.py $param
./run_any.sh ./calc_likelihood_bounds.py $param
./run_any.sh ./make_zoom_param.py $param

./run_any.sh ./core_fit $param_zoom
./run_any.sh ./plot_zmrs.py $param_zoom
./run_any.sh ./calc_likelihood_bounds.py $param_zoom
