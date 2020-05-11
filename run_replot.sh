#!/bin/bash

make_scirpt=$1
unset DISPLAY


# ./run_any_make.sh ./calc_likelihood_bounds.py params/rmba/make_all_LJ.McClintock.high_richness.low_rez.min20.sh
# ./run_any_make.sh ./calc_likelihood_bounds.py params/rmba/make_all_OR.Baxter.high_richness.low_rez.min20.sh
# ./run_any_make.sh ./calc_likelihood_bounds.py params/rmba/make_all_OR.Farahi.high_richness.low_rez.min20.sh
./run_any_make.sh ./calc_likelihood_bounds.py params/rmba/make_all_OR.McClintock.high_richness.low_rez.min20.sh
# ./run_any_make.sh ./calc_likelihood_bounds.py params/rmba/make_all_OR.high_richness.low_rez.min20.sh

# ./run_any_make.sh ./calc_likelihood_bounds.py params/rmba/make_all_LJ.McClintock.high_richness.low_rez.min20.sh zoom
# ./run_any_make.sh ./calc_likelihood_bounds.py params/rmba/make_all_OR.Baxter.high_richness.low_rez.min20.sh zoom
# ./run_any_make.sh ./calc_likelihood_bounds.py params/rmba/make_all_OR.Farahi.high_richness.low_rez.min20.sh zoom
./run_any_make.sh ./calc_likelihood_bounds.py params/rmba/make_all_OR.McClintock.high_richness.low_rez.min20.sh zoom
./run_any_make.sh ./calc_likelihood_bounds.py params/rmba/make_all_OR.high_richness.low_rez.min20.sh zoom
