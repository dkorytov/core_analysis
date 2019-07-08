#!/bin/bash


if [[ "$1" == "false" ]]
    then 
    make_param=false
    else 
    echo "we will make the params"
    echo ""
    echo "========================="
    echo ""
    make_param=true
fi

expected_comov_abundance=0.00241674353851

# Two elements for m200c and m200m
zmrh5_locs=("/home/dkorytov/phys/Ngal_sdss/data/rad_profile_mstar0_wmap7_simet_mean3/result/type1_weight1_mag1_clr1_result.hdf5" "/home/dkorytov/phys/Ngal_sdss/data/rad_profile_mstar0_wmap7_simet_crit3/result/type1_weight1_mag1_clr1_result.hdf5")
cluster_locs=("tmp_hdf5/clusters_OR_M200m@core_usage@.@radius@.hdf5" "tmp_hdf5/clusters_OR_M200c@core_usage@.@radius@hdf5")
cluster_types=("mean" "crit")

# model flavors
fit_r_mergers=("false" "false"  "true"  "true")
mi_bins_infos=("1024"    "128"   "128"    "40")
rd_bins_infos=("1"       "128"     "1"    "40")
rm_bins_infos=("1"         "1"    "32"    "20")

model_types=("mi" "rd" "rm" "rd_rm")

# fit to what
cost_abundances=("false" "true")
cost_types=("" "abund/")

# central
central_usage=("false" "true")
central_types=("" "_cen")

#Core def
core_usage=("" ".peak")
core_types=("" "_peak")

mstar_vals=("-1" "0" "0.5" "1")
for cluster_i in 0;do               # just m200m
    for mstar_i in 0 1 2 3; do      # for mstar from -1 -> +1
	for model_i in 2 3; do      # just infall and disruption
	    for cen_i in 0;  do     # no central
		for core_i in 0; do # no peak
		    param_fname="params/cfn/simet/mstar0/${cluster_types[$cluster_i]}/rscan_${radii_type[$rad_i]}.${model_types[$model_i]}${central_types[${cen_i}]}${core_types[$core_i]}.param"
		    # echo ${cluster_i} ${model_i} ${cost_i}
		    echo $param_fname
		    if [ $make_param = true ]; then
			cp params/templates/template.param $param_fname
			# SDSS zmr
			sed -i "s%@zmrh5_loc@%${zmrh5_locs[$cluster_i]}%g"         $param_fname
			# cluster loc
			sed -i "s%@cluster_loc@%${cluster_locs[$cluster_i]}%g"     $param_fname
			# cluster radius
			sed -i "s/@radius@/${radii_type[$rad_i]}/g" $param_fname
			sed -i "s/@fit_r_merger@/${fit_r_mergers[$model_i]}/g"     $param_fname
			sed -i "s/@mi_binds_info_num@/${mi_bins_infos[$model_i]}/g"  $param_fname
			sed -i "s/@rd_binds_info_num@/${rd_bins_infos[$model_i]}/g"  $param_fname
			sed -i "s/@rm_binds_info_num@/${rm_bins_infos[$model_i]}/g"  $param_fname
			sed -i "s/@cost_abundance@/false/g" $param_fname
			sed -i "s/@central_usage@/${central_usage[$cen_i]}/g" $param_fname
			sed -i "s/@core_usage@/${core_usage[$core_i]}/g" $param_fname
			sed -i "s/@cluster_radius@/${radii_type[$rad_i]}/g" $param_fname
			## global conversions not related to the for loop
			sed -i "s/@radial_bin_start@/0/g" $param_fname
			sed -i "s/@expected_comov_abundance@/$expected_comov_abundance/g" $param_fname
		    fi
		done
	    done
	done
    done
done    
