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
    bins = np.linspace(start, end, int(num))
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
            rd_factor = 1.0
            if mp == "rd":
                rd_factor = 1000.0

            old_values = base_param.get_float_list(mp+"_bins_info")
            center = model_fit_param.get_float(mp)
            lower, upper = model_fit_param.get_float_list(mp+"_limits")
            print("\n\nparam:", mp)
            old_range = (old_values[1]-old_values[0])*rd_factor
            print(" old range: {:.3f}->{:.3f}".format(old_values[0]*rd_factor, old_values[1]*rd_factor))
            print("\trange: {:.3f}".format((old_values[1]-old_values[0])*rd_factor))
            print(" fit: {:.3f} +{:.3f} -{:.3f}".format(center, upper-center, center-lower))
            print("\trange: {:3f}".format(upper-lower))

            err = np.max([center-lower, upper-center, spacing[mp]])

            # print(center-lower, upper-center)
            print(" target err: {:.3f}".format(err))
            # if err< 0.002:
            #     err = 0.002
            new_lower_limit = lower - err*3
            new_upper_limit = upper + err*3
            if new_lower_limit < 0:
                new_lower_limit = 0
            print(" New limits: {:.3f} -> {:.3f}".format(new_lower_limit, new_upper_limit))
            new_range = new_upper_limit - new_lower_limit
            print("\trange: {:.3f}".format(new_range))
            print("\tchange: {:.3f}".format(new_range/old_range))
            if mp != 'rd':
                base_param.set_var(mp+"_bins_info", "{} {} {}".format(new_lower_limit, new_upper_limit, old_values[2]))
            else:
                base_param.set_var(mp+"_bins_info", "{} {} {}".format(new_lower_limit/1000, new_upper_limit/1000, old_values[2]))

    with open(input_param_fname.replace(".param", "_zoom.param"), 'w') as f:
        f.write(str(base_param))
                   
    


if __name__ == "__main__":
    make_zoom_param(sys.argv[1])
