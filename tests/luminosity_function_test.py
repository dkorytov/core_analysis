#!/usr/bin/env python2.7


import numpy as np
import matplotlib.pyplot as plt




def schecter(phi, mstar, alpha):
    return lambda m : 0.4*np.log(10)*phi*(10**(0.4*(mstar-m)))**(alpha+1.0)*np.exp(-10.0**(0.4*(mstar-m)))


lum_funct = schecter(1.49*10**-2.0,-20.4, -1.05)


mag = np.linspace(-16,-26,100)
phi = lum_funct(mag)
plt.figure()
plt.plot(mag,phi*256**3+1.0,'x-')
plt.yscale('log')
print lum_funct(-21)

plt.show()
