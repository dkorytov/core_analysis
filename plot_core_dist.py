#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import dtk
import sys



if __name__ == "__main__":
    param = dtk.Param(sys.argv[1])
    steps = param.get_int_list("steps")
    core_loc = param.get_string("core_loc")
    for s in steps:
        print "step: ",s
        cfname = core_loc.replace("${step}",str(s))
        radius = dtk.gio_read(cfname,"radius")
        infall_mass = dtk.gio_read(cfname,"infall_mass")
        infall_step = dtk.gio_read(cfname,"infall_step")
