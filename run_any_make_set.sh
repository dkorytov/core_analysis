#!/bin/bash

make_scirpt=$1


echo $make_scirpt

./run_any_make.sh ./core_fit $make_scirpt
./run_any_make.sh ./calc_likehood_bounds.py $make_scirpt
./run_any_make.sh ./make_zoom_make_scirpt.py $make_scirpt

./run_any_make.sh ./core_fit $make_scirpt zoom
./run_any_make.sh ./calc_likehood_bounds.py $make_scirpt zoom
./run_any_make.sh ./make_zoom_make_scirpt.py $make_scirpt zoom
