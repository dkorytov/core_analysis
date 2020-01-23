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

# Two elements for m200m and m200c
zmrh5_locs=( '/data/a/cpac/dkorytov/data/Ngal_sdss_old/data/rad_profile_mstar@mstarval@_wmap7_simet_mean3/result/type1_weight1_mag1_clr1_result.hdf5'
	     '/data/a/cpac/dkorytov/data/Ngal_sdss_old/data/rad_profile_mstar@mstarval@_wmap7_simet_crit3/result/type1_weight1_mag1_clr1_result.hdf5')
#"/home/dkorytov/phys/Ngal_sdss/data/rad_profile_mstar@mstarval@_wmap7_simet_mean3/result/type1_weight1_mag1_clr1_result.hdf5" "/home/dkorytov/phys/Ngal_sdss/data/rad_profile_mstar@mstarval@_wmap7_simet_crit3/result/type1_weight1_mag1_clr1_result.hdf5")
cluster_locs=("tmp_hdf5/clusters_AQ_M200m.hdf5" "tmp_hdf5/clusters_AQ_M200c.hdf5")
cluster_types=("mean" "crit")
cluster_load_num="500"
# model flavors
fit_r_mergers=("false" "false"  "true"  "true")
mi_bins_infos=("1024"    "128"   "128"    "32")
rd_bins_infos=("1"       "128"     "1"    "32")
rm_bins_infos=("1"         "1"    "32"    "16")

model_types=("mi" "rd" "rm" "rd_rm")

# fit to what
cost_abundances=("false" "true")
cost_types=("" "abund/")

# mstars=("-1" "0" "0.5" "1")
mstars=("0")
displacements=("0.00" "0.01" "0.02" "0.05" "0.10" "0.20" "0.50" "1.00")
for mstar in "${mstars[@]}"; do
    for cluster_i in 1; do
	for model_i in 1 ; do 
	    for cost_i in 0; do
		for displacement in "${displacements[@]}"; do
		    prefix="AQ_${displacement}"
		    param_fname="params/rmba/simet/${cluster_types[$cluster_i]}/mstar${mstar}/${cost_types[$cost_i]}/center_displace/${prefix}_${model_types[$model_i]}.param"
		    dir=${param_fname%/*}
		    # echo ${cluster_i} ${model_i} ${cost_i}
		    echo $param_fname
		    if [ $make_param = true ]; then
			mkdir -p $dir
			# if we are mstar -1 params or it's mstar0 and mi
			# model use template_mstar-1.param template file
			# for the expanded M_infall range
			if [ $mstar == -1 ] || { [ $model_i == 0 ] && [ $mstar == 0 ]; }
			then
			    cp params/templates/template_mstar-1.param $param_fname
			else 
			    cp params/templates/template.param $param_fname
			fi
			sed -i "s%@zmrh5_loc@%${zmrh5_locs[$cluster_i]}%g"         $param_fname
			sed -i "s/@cluster_load_num@/${cluster_load_num}/g"        $param_fname   
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
			sed -i "s/@radial_bin_start@/0/g" $param_fname

			## Displacement parameters
			sed -i "s/displace_cluster_center_distance 0.0/displace_cluster_center_distance ${displacement}/g" $param_fname
			sed -i "s/force_center_on_central false/force_center_on_central true/g" $param_fname

		    fi
		done
	    done
	done
    done
done
