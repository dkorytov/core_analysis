#!/bin/bash

make_scirpt=$1
unset DISPLAY
echo "Running set with 1 level zoom"
echo $make_scirpt
echo ""
date
echo ""

./run_any_make.sh ./core_fit $make_scirpt
./run_any_make.sh ./calc_likelihood_bounds.py $make_scirpt
./run_any_make.sh ./make_zoom_param.py $make_scirpt
echo ""
echo "Phase 1 Done:"
date
echo ""

./run_any_make.sh ./core_fit $make_scirpt zoom
./run_any_make.sh ./calc_likelihood_bounds.py $make_scirpt zoom
./run_any_make.sh ./make_zoom_param.py $make_scirpt zoom
echo ""
echo "Phase 2 Done:"
date
echo ""

./run_any_make.sh ./core_fit $make_scirpt zoom zoom
./run_any_make.sh ./calc_likelihood_bounds.py $make_scirpt zoom zoom 
./run_any_make.sh ./make_zoom_param.py $make_scirpt zoom zoom

echo ""
echo "All done"
date
echo ""
