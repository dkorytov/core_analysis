#!/bin/bash


if [[ "$1" == "false" ]]; then 
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
cluster_locs=("tmp_hdf5/clusters_OR_M200m.peak.hdf5" "tmp_hdf5/clusters_OR_M200c.peak.hdf5")
cluster_types=("mean" "crit")

# model flavors
fit_r_mergers=("false" "false"  "true"  "true")
mi_bins_infos=("256"   "128"    "128"   "40"  )
rd_bins_infos=("1"     "128"      "1"   "40"  )
rm_bins_infos=("1"       "1"     "32"   "20"  )

model_types=("mi" "rd" "rm" "rd_rm")

# fit to what
cost_abundances=("false" "true")
cost_types=("" "abund/")

for cluster_i in 0 1;do
    for model_i in 0 1 2 3; do 
	for cost_i in 0 1; do
	    param_fname="params/cfn/simet/mstar0/${cluster_types[$cluster_i]}/${cost_types[$cost_i]}a_${model_types[$model_i]}_peak.param"
	    # echo ${cluster_i} ${model_i} ${cost_i}
	    echo $param_fname
	    if [ $make_param = true ]; then
		cp params/templates/template.param $param_fname
		sed -i "s%@zmrh5_loc@%${zmrh5_locs[$cluster_i]}%g"         $param_fname
		sed -i "s%@cluster_loc@%${cluster_locs[$cluster_i]}%g"     $param_fname
		sed -i "s/@fit_r_merger@/${fit_r_mergers[$model_i]}/g"     $param_fname
		sed -i "s/@mi_binds_info_num@/${mi_bins_infos[$model_i]}/g"  $param_fname
		sed -i "s/@rd_binds_info_num@/${rd_bins_infos[$model_i]}/g"  $param_fname
		sed -i "s/@rm_binds_info_num@/${rm_bins_infos[$model_i]}/g"  $param_fname
		sed -i "s/@cost_abundance@/${cost_abundances[$cost_i]}/g" $param_fname
		## global conversions not related to the for loop
		sed -i "s/@expected_comov_abundance@/${expected_comov_abundance}/g" $param_fname
		## global conversions not related to the for loop
		sed -i "s/@radial_bin_start@/0/g" $param_fname
	    done
	done
    done
done
