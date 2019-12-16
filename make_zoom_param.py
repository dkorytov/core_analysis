#!/usr/bin/env python3


from __future__ import print_function, division 
import numpy as np
import os




import dtk
import sys

def get_spacing(bin_info):
    start = bin_info[0]
    end = bin_info[1]
    num = bin_info[2]
    if num < 2:
        return 0
    bins = np.linspace(start, end, num)
    return np.diff(bins)[0]

def make_zoom_param(input_param_fname):
    base_param = dtk.Param(input_param_fname)
    mi_bins_info = base_param.get_float_list("mi_bins_info")
    rd_bins_info = base_param.get_float_list("rd_bins_info")
    rm_bins_info = base_param.get_float_list("rm_bins_info")
    spacing = {}
    spacing["mi"] = get_spacing(mi_bins_info)
    spacing["rd"] = get_spacing(rd_bins_info)
    spacing["rm"] = get_spacing(rm_bins_info)
    model_fit_param = dtk.Param("figs/"+input_param_fname+"/calc_likelihood_bounds.py/grid_fit_param.txt")
    model_params = ['mi', 'rd', 'rm']
    for mp in model_params:
        if mp in model_fit_param:
            center = model_fit_param.get_float(mp)
            lower, upper = model_fit_param.get_float_list(mp+"_limits")
            print("center: ", center)
            print("limits: ", lower, upper)
            old_values = base_param.get_float_list(mp+"_bins_info")
            err = np.max([center-lower, upper-center, spacing[mp]])
            print(mp)
            print(center-lower, upper-center)
            print(err)
            # if err< 0.002:
            #     err = 0.002
            new_lower_limit = lower - err*3
            new_upper_limit = upper + err*3
            if new_lower_limit < 0:
                new_lower_limit = 0
            print("New limits: ", new_lower_limit, new_upper_limit)
            if mp != 'rd':
                base_param.set_var(mp+"_bins_info", "{} {} {}".format(new_lower_limit, new_upper_limit, old_values[2]))
            else:
                base_param.set_var(mp+"_bins_info", "{} {} {}".format(new_lower_limit/1000, new_upper_limit/1000, old_values[2]))

    with open(input_param_fname.replace(".param", "_zoom.param"), 'w') as f:
        f.write(str(base_param))
                   
    
    print(base_param.__str__())

if __name__ == "__main__":
    make_zoom_param(sys.argv[1])
