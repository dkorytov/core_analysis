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
zmrh5_locs=("/home/dkorytov/phys/Ngal_sdss/data/rad_profile_mstar@mstarval@_wmap7_spider_mean/result/type1_weight1_mag1_clr1_result.hdf5" "/home/dkorytov/phys/Ngal_sdss/data/rad_profile_mstar@mstarval@_wmap7_spider/result/type1_weight1_mag1_clr1_result.hdf5")
cluster_locs=("tmp_hdf5/clusters_OR_M200m.hdf5" "tmp_hdf5/clusters_OR_M200c.hdf5")
cluster_types=("mean" "crit")

# model flavors
fit_r_mergers=("false" "false"  "true"  "true")
mi_bins_infos=("1024"    "128"   "128"    "32")
rd_bins_infos=("1"       "128"     "1"    "32")
rm_bins_infos=("1"         "1"    "32"    "16")

model_types=("mi" "rd" "rm" "rd_rm")

# fit to what
cost_abundances=("false" "true")
cost_types=("" "abund/")

mstars=("-1" "0" "0.5" "1")

for mstar in "${mstars[@]}"; do
    for cluster_i in 0 ;do
	for model_i in 0 1 ; do 
	    for cost_i in 0; do
		param_fname="params/cfn/simet/mstar${mstar}/${cluster_types[$cluster_i]}/${cost_types[$cost_i]}spider_${model_types[$model_i]}.param"
		# echo ${cluster_i} ${model_i} ${cost_i}
		echo $param_fname
		if [ $make_param = true ]; then
		    # if we are mstar -1 params or it's mstar0 and mi
		    # model use template_mstar-1.param template file
		    # for the expanded M_infall range
		    if [ $mstar == -1 ] || { [ $model_i == 0 ] && [ $mstar == 0 ]; } then
			cp params/templates/template_mstar-1.param $param_fname
		    else 
			cp params/templates/template.param $param_fname
		    fi
		    sed -i "s%@zmrh5_loc@%${zmrh5_locs[$cluster_i]}%g"         $param_fname
		    sed -i "s%@mstarval@%${mstar}%g"                           $param_fname
		    sed -i "s%@cluster_loc@%${cluster_locs[$cluster_i]}%g"     $param_fname
		    sed -i "s/@fit_r_merger@/${fit_r_mergers[$model_i]}/g"     $param_fname
		    sed -i "s/@mi_binds_info_num@/${mi_bins_infos[$model_i]}/g"  $param_fname
		    sed -i "s/@rd_binds_info_num@/${rd_bins_infos[$model_i]}/g"  $param_fname
		    sed -i "s/@rm_binds_info_num@/${rm_bins_infos[$model_i]}/g"  $param_fname
		    sed -i "s/@cost_abundance@/${cost_abundances[$cost_i]}/g" $param_fname
		    ## global conversions not related to the for loop
		    sed -i "s/@expected_comov_abundance@/${expected_comov_abundance}/g" $param_fname
		    ## global conversions not related to the for loop
		    sed -i "s/@radial_bin_start@/2/g" $param_fname
		    sed -i "s/cost_min_cluster_num  0/cost_min_cluster_num  15/g" $param_fname
		fi
	    done
	done
    done
done
