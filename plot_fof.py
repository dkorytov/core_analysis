#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import dtk



param = dtk.Param("tmp/fof_test.param")
x = np.array(param.get_float_list("x"))
y = np.array(param.get_float_list("y"))
z = np.array(param.get_float_list("z"))
c = np.array(param.get_int_list("c"))
clrs =  np.unique(c)
print clrs.size,'/',c.size

plt.figure()
plt.scatter(x,y,c=c,edgecolor='none')
plt.xlim((0,5))
plt.ylim((0,5))
plt.grid()
plt.show()
