#!/usr/bin/env python2.7

from __future__ import print_function, division

from matplotlib import rc
rc('font',**{'family': 'serif',
             'serif':  ['DejaVu'],
             'size':   15})
rc('text', usetex=True)

import dtk
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import numpy as np




def plot(val_list):
    plt.figure()
    y = np.arange(len(val_list))[::-1]
    name_list = []
    for i, (name, vals) in enumerate(val_list):
        print(i, name, vals)
        name_list.append(name)
        print(vals[0], y[i], vals[0]-vals[1], vals[2]-vals[1])
        plt.plot(vals[0], y[i], 'x')
        plt.errorbar(vals[0], y[i], fmt='o', xerr = [[vals[0]-vals[1]], [vals[2]-vals[1]]])
    plt.ylim((np.min(y)-0.5, np.max(y)+0.5))
    # plt.xlim([12.1, 12.3])
    plt.show()


# data, model, parameter, values

def load_fit_limits(fname):
    if fname in load_fit_limits.cache:
        return load_fit_limits.cache[fname]
    pfile = dtk.Param(fname)
    result = {}
    for p in ['mi', 'rd', 'rm']:
        if p in pfile:
            if p == 'rd': #rescale rd into kpc
                f =1000
            else:
                f =1
            result[p] = pfile.get_float(p)*f
            result[p+'_limits'] = np.array(pfile.get_float_list(p+'_limits'))*f
            result[p+'_lower_err'] = result[p] - result[p+'_limits'][0]
            result[p+'_upper_err'] = result[p+'_limits'][1] - result[p]
    if 'X_red' in pfile:
        result['x2'] = pfile.get_float('X_red')
        print(result['x2'])
    else:
        result['x2'] = np.nan
    return result

load_fit_limits.cache = {}

def load_all():
    pattern = "figs/params/cfn/simet/{}{}/{}.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    models = {"rd3": "D", "rm3": "M", "rd_rm3": "M + D"}
    mass_defs = {"crit":"M200c", "mean":"M200m" }
    fit_types = {"":"prof only",  "/abund": "prof + abund"}
    data_num =  len(mass_defs)*len(fit_types)
    data_clr = ['b', 'g', 'r', 'c']
    plt.figure()
    for model_i, model in enumerate(models.keys()):
        offset=.25
        clr_i = 0
        for mass_def in mass_defs.keys():
            for fit_type in fit_types.keys():
                if(mass_def == "crit"  and fit_type == "/abund"):
                    print("hhhmmm")
                    continue
                print(model, mass_def, fit_type)
                res = load_fit_limits(pattern.format(mass_def, fit_type, model))
                plt.plot(res['mi'], model_i+offset, "o", color=data_clr[clr_i])
                plt.errorbar(res['mi'], model_i+offset, xerr=[[res['mi_lower_err']], [res['mi_upper_err']]], color=data_clr[clr_i])
                offset+=0.5/float(data_num-1)
                clr_i +=1

def plot_mstar_minus1():
    pattern = "figs/params/cfn/simet/mstar-1/{}{}/{}_zoom.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    models = ["a_rd", "a_rm", "a_rd_rm"]
    data_inputs = [['mean', ''],
                   ['mean', '/abund'],
                   ['crit', ''],
                   ['crit', '/abund']]
    data_input_labels = ['M200m, profile', 
                         'M200m, profile + abund',
                         'M200c, profile',
                         'M200c, profile + abund']
    data_clr = ['b', 'b', 'r', 'r']
    data_mfc = ['b', 'none', 'r','none']
    load_all2(pattern, models, data_inputs, data_input_labels, data_clr, data_mfc, title="Mstar-1")

def plot_mstar0():
    pattern = "figs/params/cfn/simet/{}{}/{}_zoom.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    models = ["a_rd", "a_rm", "a_rd_rm"]
    data_inputs = [['mean', ''],
                   ['mean', '/abund'],
                   ['crit', ''],
                   ['crit', '/abund']]
    data_input_labels = ['M200m, profile', 
                         'M200m, profile + abund',
                         'M200c, profile',
                         'M200c, profile + abund']
    data_clr = ['b', 'b', 'r', 'r']
    data_mfc = ['b', 'none', 'r','none']
    load_all2(pattern, models, data_inputs, data_input_labels, data_clr, data_mfc, title="Mstar+0")

def plot_mstar1():
    pattern = "figs/params/cfn/simet/mstar1/{}{}/{}_zoom.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    models = ["a_rd", "a_rm", "a_rd_rm"]
    data_inputs = [['mean', ''],
                   ['mean', '/abund'],
                   ['crit', ''],
                   ['crit', '/abund']]
    data_input_labels = ['M200m, profile', 
                         'M200m, profile + abund',
                         'M200c, profile',
                         'M200c, profile + abund']
    data_clr = ['b', 'b', 'r', 'r']
    data_mfc = ['b', 'none', 'r','none']
    load_all2(pattern, models, data_inputs, data_input_labels, data_clr, data_mfc, title="Mstar+1")

def plot_mstar05():
    pattern = "figs/params/cfn/simet/mstar0.5/{}{}/{}_zoom.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    models = ["a_rd", "a_rm", "a_rd_rm"]
    data_inputs = [['mean', ''],
                   ['mean', '/abund'],
                   ['crit', ''],
                   ['crit', '/abund']]
    data_input_labels = ['M200m, profile', 
                         'M200m, profile + abund',
                         'M200c, profile',
                         'M200c, profile + abund']
    data_clr = ['b', 'b', 'r', 'r']
    data_mfc = ['b', 'none', 'r','none']
    load_all2(pattern, models, data_inputs, data_input_labels, data_clr, data_mfc, title="Mstar+0.5")

def plot_mstar05b():
    pattern = "figs/params/cfn/simet/mstar0.5/{}{}/{}_zoom.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    models = ["ab_rd", "ab_rm", "ab_rd_rm"]
    data_inputs = [['mean', ''],
                   ['mean', '/abund'],
                   ['crit', ''],
                   ['crit', '/abund']]
    data_input_labels = ['M200m, profile', 
                         'M200m, profile + abund',
                         'M200c, profile',
                         'M200c, profile + abund']
    data_clr = ['b', 'b', 'r', 'r']
    data_mfc = ['b', 'none', 'r','none']
    load_all2(pattern, models, data_inputs, data_input_labels, data_clr, data_mfc, title="Mstar+0.5 (b)")

def plot_mstar0_pc():
    file_patern = "figs/params/cfn/simet/mstar0/${mass_def}/pc_${model}${cent}${peak}.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    mass_defs = {'mean', 'crit'}
    model_defs = {'rd', 'rm', 'rd_rm'}
    cent_defs = {'', '_cen'}
    peak_defs = {'', '_peak'}
    load_pc(file_patern, mass_defs, model_defs, cent_defs, peak_defs)

def plot_mstar05_pc():
    file_patern = "figs/params/cfn/simet/mstar0.5/${mass_def}/pc_${model}${cent}${peak}.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    mass_defs = {'mean', 'crit'}
    model_defs = {'rd', 'rm', 'rd_rm'}
    cent_defs = {'', '_cen'}
    peak_defs = {'', '_peak'}
    load_pc(file_patern, mass_defs, model_defs, cent_defs, peak_defs)

def plot_mstar1_pc():
    file_patern = "figs/params/cfn/simet/mstar1/${mass_def}/pc_${model}${cent}${peak}.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    mass_defs = {'mean', 'crit'}
    model_defs = {'rd', 'rm', 'rd_rm'}
    cent_defs = {'', '_cen'}
    peak_defs = {'', '_peak'}
    load_pc(file_patern, mass_defs, model_defs, cent_defs, peak_defs)

def plot_mstar0_QC():
    file_patern = "figs/params/cfn/simet/mstar0/${mass_def}/qc_${model}${cent}${peak}.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    mass_defs = {'mean', 'crit'}
    model_defs = {'rd', 'rm', 'rd_rm'}
    cent_defs = {'',}
    peak_defs = {'',}
    load_pc(file_patern, mass_defs, model_defs, cent_defs, peak_defs)

def plot_mstar05_QC():
    file_patern = "figs/params/cfn/simet/mstar0.5/${mass_def}/qc_${model}${cent}${peak}.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    mass_defs = {'mean', 'crit'}
    model_defs = {'rd', 'rm', 'rd_rm'}
    cent_defs = {'',}
    peak_defs = {'',}
    load_pc(file_patern, mass_defs, model_defs, cent_defs, peak_defs)

def plot_mstar1_QC():
    file_patern = "figs/params/cfn/simet/mstar1/${mass_def}/qc_${model}${cent}${peak}.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    mass_defs = {'mean', 'crit'}
    model_defs = {'rd', 'rm', 'rd_rm'}
    cent_defs = {'',}
    peak_defs = {'',}
    load_pc(file_patern, mass_defs, model_defs, cent_defs, peak_defs)

def plot_mstar0_OR():
    pattern = "figs/params/cfn/simet/mstar0/{}{}/{}.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    models = ["a2_mi", "a2_rd", "a2_rm", "a2_rd_rm"]
    data_inputs = [['mean', '']]
    data_input_labels = ['M200m, profile']
    data_clr = ['b', ]
    data_mfc = ['b', ]
    load_all3(pattern, models, data_inputs, data_input_labels, data_clr, data_mfc, title=None)


def load_all2(pattern, models, data_inputs, data_input_labels, data_clr, data_mfc, title=None):
    plt.figure()
    if title is not None:
        plt.suptitle(title)
    gs = gridspec.GridSpec(5,4, hspace=0.1)
    model_params_ax = {'mi':plt.subplot(gs[1:,0]),
                       'rd':plt.subplot(gs[1:,1]),
                       'rm':plt.subplot(gs[1:,2]),
                       'x2':plt.subplot(gs[1:,3])}
    model_params = ['mi', 'rd', 'rm']
    for model_i, model in enumerate(models):
        for data_i in range(0, len(data_inputs)):
            model_param_fit = load_fit_limits(pattern.format(data_inputs[data_i][0], data_inputs[data_i][1], model))
            print("\n\n", model, data_inputs[data_i])
            if len(data_inputs) == 1:
                y = -model_i
            else:
                y=-( model_i -0.15 + 0.3/(len(data_inputs)-1)*data_i)
            for model_param in model_params:
                if model_param in model_param_fit:
                    print(model_param)
                    ax = model_params_ax[model_param]
                    x = model_param_fit[model_param]
                    xerr = [[model_param_fit[model_param+"_lower_err"]], [model_param_fit[model_param+"_upper_err"]]]
                    ax.errorbar(x, y, xerr=xerr, fmt='o', color=data_clr[data_i], mfc=data_mfc[data_i], mec = data_clr[data_i])
            model_params_ax['x2'].errorbar(model_param_fit['x2'], y, fmt='o', color=data_clr[data_i], mfc=data_mfc[data_i], mec=data_clr[data_i])
    model_params_ax['rd'].set_ylim(model_params_ax['mi'].set_ylim())
    model_params_ax['rm'].set_ylim(model_params_ax['mi'].set_ylim())
    model_params_ax['x2'].set_ylim(model_params_ax['mi'].set_ylim())
    model_params_ax['mi'].set_xlabel('$M_{infall}$')
    model_params_ax['rd'].set_xlabel('$R_{disrupt}$')
    model_params_ax['rm'].set_xlabel('$M_{merge}$')
    model_params_ax['x2'].set_xlabel(r'$\chi^{2}_{reduced}$')
    for (key, ax) in model_params_ax.iteritems():
        # if key != 'x2':
        #     ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
        # else:
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax.yaxis.set_ticks([-0,-1,-2, -3])
        ax.yaxis.set_ticklabels(["Only\nDisruption", "Only\nMerging", "Disruption\n+\nMerging"])
        for tick in ax.xaxis.get_ticklabels():
            tick.set_rotation(45)


    model_params_ax['rd'].yaxis.set_ticklabels([])
    model_params_ax['rm'].yaxis.set_ticklabels([])
    model_params_ax['x2'].yaxis.set_ticklabels([])
    # model_params_ax['rd'].xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax = model_params_ax['rd']
    ax.plot([],[], linewidth=5, color='b', label='$M_{200m}$')
    ax.plot([],[], linewidth=5, color='r', label='$M_{200c}$')
    ax.errorbar([],[],xerr=[], fmt='o', color='k', mec='k', mfc='k',
                label='profile only')
    ax.errorbar([],[],xerr=[], fmt='o', color='k', mec='k', mfc='none',
                label='profile + abundance')
    ax.legend(bbox_to_anchor=(-1., 1.02, 4, .102), loc='lower center',
           ncol=2, mode="expand", borderaxespad=0.)
    
    dtk.save_figs("figs/"+__file__+"/")

def load_all3(pattern, models, data_inputs, data_input_labels, data_clr, data_mfc, title=None):
    plt.figure(figsize=(5,3))
    if title is not None:
        plt.suptitle(title)
    gs = gridspec.GridSpec(1,4, hspace=0.1)
    model_params_ax = {'mi':plt.subplot(gs[:,0]),
                       'rd':plt.subplot(gs[:,1]),
                       'rm':plt.subplot(gs[:,2]),
                       'x2':plt.subplot(gs[:,3])}
    model_params = ['mi', 'rd', 'rm']
    for model_i, model in enumerate(models):
        for data_i in range(0, len(data_inputs)):
            model_param_fit = load_fit_limits(pattern.format(data_inputs[data_i][0], data_inputs[data_i][1], model))
            print("\n\n", model, data_inputs[data_i])
            if len(data_inputs) == 1:
                y = -model_i
            else:
                y=-( model_i -0.15 + 0.3/(len(data_inputs)-1)*data_i)
            for model_param in model_params:
                if model_param in model_param_fit:
                    print(model_param)
                    ax = model_params_ax[model_param]
                    x = model_param_fit[model_param]
                    xerr = [[model_param_fit[model_param+"_lower_err"]], [model_param_fit[model_param+"_upper_err"]]]
                    ax.errorbar(x, y, xerr=xerr, fmt='o', color=data_clr[data_i], mfc=data_mfc[data_i], mec = data_clr[data_i])
            # model_params_ax['x2'].errobar(model_param_fit['x2'], y, fmt='o', color=data_clr[data_i], mfc=data_mfc[data_i], mec=data_clr[data_i])
            model_params_ax['x2'].barh(y, model_param_fit['x2'], align='center', height=0.2, edgecolor='none' )
    model_params_ax['mi'].set_ylim([ -3.2, 0.2,])
    model_params_ax['rd'].set_ylim(model_params_ax['mi'].set_ylim())
    model_params_ax['rm'].set_ylim(model_params_ax['mi'].set_ylim())
    model_params_ax['x2'].set_ylim(model_params_ax['mi'].set_ylim())
    model_params_ax['mi'].set_xlabel('$M_{infall}$')
    model_params_ax['rd'].set_xlabel('$R_{disrupt}$')
    model_params_ax['rm'].set_xlabel('$M_{merge}$')
    model_params_ax['x2'].set_xlabel(r'$\chi^{2}_{reduced}$')
    for (key, ax) in model_params_ax.iteritems():
        # if key != 'x2':
        #     ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
        # else:
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax.yaxis.set_ticks([-0,-1,-2, -3])
        ax.yaxis.set_ticklabels(["Mi", "Rd", "Rm", "RdRm"])
        for tick in ax.xaxis.get_ticklabels():
            tick.set_rotation(45)


    model_params_ax['rd'].yaxis.set_ticklabels([])
    model_params_ax['rm'].yaxis.set_ticklabels([])
    model_params_ax['x2'].yaxis.set_ticklabels([])
    # model_params_ax['rd'].xaxis.set_major_locator(ticker.MultipleLocator(5))
    # ax = model_params_ax['rd']
    # ax.plot([],[], linewidth=5, color='b', label='$M_{200m}$')
    # ax.plot([],[], linewidth=5, color='r', label='$M_{200c}$')
    # ax.errorbar([],[],xerr=[], fmt='o', color='k', mec='k', mfc='k',
    #             label='profile only')
    # ax.errorbar([],[],xerr=[], fmt='o', color='k', mec='k', mfc='none',
    #             label='profile + abundance')
    # ax.legend(bbox_to_anchor=(-1., 1.02, 4, .102), loc='lower center',
    #        ncol=2, mode="expand", borderaxespad=0.)
    # plt.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.3,left=0.15)


    dtk.save_figs("figs/"+__file__+"/", extension='.pdf')

def load_pc(pattern, mass_defs, model_defs, cent_defs, peak_defs, title=None):
    plt.figure()
    if title is not None:
        plt.suptitle(title)
    gs = gridspec.GridSpec(5,4, hspace=0.1)
    model_params_ax = {'mi':plt.subplot(gs[1:,0]),
                       'rd':plt.subplot(gs[1:,1]),
                       'rm':plt.subplot(gs[1:,2]),
                       'x2':plt.subplot(gs[1:,3])}
    model_params = ['mi', 'rd', 'rm']
    index_max = len(mass_defs)*len(cent_defs)*len(peak_defs)-1
    for mass_i, mass in enumerate(mass_defs):
        data_mfc = [None, 'none'][mass_i]
        for model_i, model in enumerate(model_defs):
            for cent_i, cent in enumerate(cent_defs):
                for peak_i, peak in enumerate(peak_defs):
                    fname = pattern.replace("${mass_def}", mass).replace("${model}", model).replace("${cent}", cent).replace("${peak}", peak)
                    print(fname)
                    model_param_fit = load_fit_limits(fname)
                    index = peak_i + cent_i*len(peak_defs) + mass_i*len(cent_defs)*len(peak_defs)
                    y = 3-model_i - 0.3*(index/index_max-0.5)
                    color = ['b', 'r', 'g', 'm'][cent_i + 2*peak_i]
                    for model_param in model_params:
                        if model_param in model_param_fit:
                            print(model_param)
                            ax = model_params_ax[model_param]
                            x = model_param_fit[model_param]
                            xerr = [[model_param_fit[model_param+"_lower_err"]], [model_param_fit[model_param+"_upper_err"]]]
                            ax.errorbar(x, y, xerr=xerr, fmt='o', mfc=data_mfc, mec = color, color = color)
                    model_params_ax['x2'].errorbar(model_param_fit['x2'], y, fmt='o',mfc=data_mfc, mec = color, color=color)
    model_params_ax['rd'].set_ylim(model_params_ax['mi'].set_ylim())
    model_params_ax['rm'].set_ylim(model_params_ax['mi'].set_ylim())
    model_params_ax['x2'].set_ylim(model_params_ax['mi'].set_ylim())
    model_params_ax['mi'].set_xlabel('$M_{infall}$')
    model_params_ax['rd'].set_xlabel('$R_{disrupt}$')
    model_params_ax['rm'].set_xlabel('$M_{merge}$')
    model_params_ax['x2'].set_xlabel(r'$\chi^{2}_{reduced}$')
    for (key, ax) in model_params_ax.iteritems():
        # if key != 'x2':
        #     ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
        # else:
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax.yaxis.set_ticks([3,2,1])
        ax.yaxis.set_ticklabels(["Only\nDisruption", "Only\nMerging", "Disruption\n+\nMerging"])
        for tick in ax.xaxis.get_ticklabels():
            tick.set_rotation(45)


    model_params_ax['rd'].yaxis.set_ticklabels([])
    model_params_ax['rm'].yaxis.set_ticklabels([])
    model_params_ax['x2'].yaxis.set_ticklabels([])
    # model_params_ax['rd'].xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax = model_params_ax['rd']
    ax.errorbar([], [], fmt='o', color='k', mfc='k', label='M200m')
    ax.errorbar([], [], fmt='o', color='k', mfc='none', label='M200c')
    ax.plot([],[], 'ob', label='Default')
    ax.plot([],[], 'or', label='central')
    ax.plot([],[], 'og', label='peak')
    ax.plot([],[], 'om', label='central+peak')
    ax.legend(bbox_to_anchor=(-1., 1.02, 4, .102), loc='lower center',
           ncol=2, mode="expand", borderaxespad=0.)
    


    

    
if __name__ == "__main__":
    # load_all()
    # plot_mstar_minus1()
    # plot_mstar0()
    # plot_mstar05()
    # plot_mstar1()
    # plot_mstar0_pc()
    # plot_mstar05_pc()
    # plot_mstar1_pc()
    # plot_mstar0_QC()
    # plot_mstar05_QC()
    # plot_mstar1_QC()
    plot_mstar0_OR()
    plt.show()
    exit()
    a = [("D", [12.2, 12.1, 12.3]),
         ("M", [12.25, 12.12, 13.0]),
         ("D+M",[12.22, 12.17, 12.25])]
    plot(a)
