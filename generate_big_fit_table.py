#!/usr/bin/env python3
"""This module takes in a "make param script", and generates a latex
table that has the best fit values. Core_fit and
calc_likelihood_bounds.py must have been run before hand. The modules
automatically attempts to read in the most "zoomed" parameter.

"""
import sys
import subprocess
from os import path

# import numpy as np
def get_zoomed_param(param_fname):
    while(path.exists(param_fname.replace(".param", "_zoom.param"))):
        param_fname = param_fname.replace(".param", "_zoom.param")
    return param_fname
def get_mstar(param_fname):
    sections = param_fname.split("/")
    mstar_substring = sections[4]
    if "mstar-1" in mstar_substring:
        return -1.0
    elif "mstar-0.5" in mstar_substring:
        return -0.5
    elif "mstar0.5" in mstar_substring:
        return 0.5
    elif "mstar0" in mstar_substring:
        return 0.0
    elif "mstar1" in mstar_substring:
        return 1.0
    else:
        raise

def get_model_flavor(param_fname):
    sections = param_fname.split("/")
    model_section = sections[-1]
    if "mi" in model_section:
        return "mi"
    elif "rd_rm" in model_section:
        return "rd_rm"
    elif "rd" in model_section:
        return "rd"
    elif "rm" in model_section:
        return "rm"
    else:
        print(param_fname)
        raise

def get_fit_file(param_fname):
    return "figs/{}/calc_likelihood_bounds.py/grid_fit_param.txt".format(param_fname)

def extract_fits_from_file(fit_fname):
    result = {}
    f = open(fit_fname)
    lines = f.readlines()
    for index in range(len(lines)//2):
        i = index*2
        if 'mi' in lines[i]:
            result['mi'] = [float(lines[i].split()[1]),
                            float(lines[i+1].split()[1]),
                            float(lines[i+1].split()[2])]
        elif 'rd' in lines[i]:
            result['rd'] = [float(lines[i].split()[1]),
                            float(lines[i+1].split()[1]),
                            float(lines[i+1].split()[2])]
        elif 'rm' in lines[i]:
            result['rm'] = [float(lines[i].split()[1]),
                            float(lines[i+1].split()[1]),
                            float(lines[i+1].split()[2])]
        elif 'cost' in lines[i]:
            result['x2'] = float(lines[i+1].split()[1])
    return result

def get_all_fits(make_param_script):
    mstar_dict ={}
    model_dict ={}
    result = subprocess.run(['bash', make_param_script, 'false'], stdout=subprocess.PIPE)
    params = result.stdout.decode("utf-8").splitlines()
    for param in params:
        print(param)
        param = get_zoomed_param(param)
        mstar = get_mstar(param)
        model_flavor = get_model_flavor(param)
        print(mstar, model_flavor)
        fit = extract_fits_from_file(get_fit_file(param))
        if mstar not in mstar_dict:
            mstar_dict[mstar] = {}
        mstar_dict[mstar][model_flavor] = fit
        if model_flavor not in model_dict:
            model_dict[model_flavor] = {}
        model_dict[model_flavor][mstar] = fit
        # print(fit)
        print('')
    return mstar_dict, model_dict

def mstar_to_lstar(val):
    return 2.5**val

def make_mstar_string(mstars):
    result = ""
    for mstar in mstars:
        result += "{:.2f} & ".format(mstar_to_lstar(mstar))
    return result[0:-2]

def get_model_table_string(model_dict, model):
    result = "\\hline\n"
    added_model = False
    parameters = ['mi', 'rd', 'rm']
    parameter_names = {'mi':'M$_{infall}$', 'rd':'R$_{disrupt}$', 'rm':'R$_{merge}$'}
    model_names = {'mi': 'Mi', 'rd':'Rd', 'rm':'Rm', 'rd_rm':'RdRm'}
    model_param_num ={'mi':1, 'rd':2, 'rm':2, 'rd_rm':3}
    for param in parameters:
        if param not in model_dict[model][0]:
            continue
        if not added_model:
            result += "\multirow{{{}}}{{*}}{{{}}} &".format(model_param_num[model], model_names[model])
            added_model = True
        else:
            result += " & "
        result += parameter_names[param] 
        for mstar in sorted(model_dict[model])[::-1]:
            vals = model_dict[model][mstar][param]
            best_fit = vals[0]
            upper_err = vals[2]-vals[0]
            lower_err = vals[0]-vals[1]
            if param != 'rd':
                result += " & ${:.3f}^{{{:+.3f}}}_{{-{:.3f}}}$".format(best_fit, upper_err, lower_err)
            else:
                result += " & ${:.1f}^{{{:+.1f}}}_{{-{:.1f}}}$".format(best_fit, upper_err, lower_err)

        result += " \\\\\n"
    return result

def generate_big_fit_table(make_param_script):
    mstar_fits, model_fits = get_all_fits(make_param_script)
    print(mstar_fits.keys())
    print(model_fits.keys())
    mstars = mstar_fits.keys()
    mstar_string = make_mstar_string(mstars)

    result = ""

    result += "\\begin{{tabular}}{{cc{}}}\n".format("c"*len(mstars))
    result += " & & \multicolumn{{{}}}{{c}}{{Galaxy Luminosity}} \\\\\n".format(len(mstars))
    result += "\cline{3-7} \n"
    result += "model & parameter & {}\\\\\n\\hline\n".format(mstar_string)
    for model in model_fits: 
        result += get_model_table_string(model_fits, model)
    result += "\\end{tabular}\n"
    print(result)

def plot_fits(make_param_script):
    mastar_fits, model_fits = get_all_fits(make_param_script)
    
if __name__ == "__main__":
    make_param_script = sys.argv[1]
    generate_big_fit_table(make_param_script)
