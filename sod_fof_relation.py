
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import cosmolopy.perturbation as prtrb 
import dtk




##Defined Parameters
#....yes,yes-I know these should be in a file

del_iso_old = 22.15
b=0.168
Hubble=0.71
Omega_DM = 0.22
Omega_BM = 0.02258/Hubble**2
Omega_M  = Omega_DM+Omega_BM
Omega_L = 1-Omega_M
prtcl_mass = 1491.**3/3200.**3*(Omega_M*2.77536627e11)

def set_prtcl_mass(mass):
    global prtcl_mass
    prtcl_mass = mass

def set_b(b_new):
    global b
    b = b_new





print "%.4g" % prtcl_mass

def del_so(del_0, a):
    bbb = del_0*(Omega_M+Omega_L*a**3)
    return bbb

def del_iso(b,a):
    return 3./(2.*np.pi*b**3)*Omega_M/a**3

def get_A(c):
    return np.log(1.0+c)-c/(1.0+c) #eq 8

def conc_old(m200,a):
    return 5.71*(m200/2e12)**-0.084*(1/a)**-0.47
 
def conc_from_v(v,z):
    return prtrb.fgrowth(z,Omega_M)**0.54*5.9*v**-0.35

def conc(m200, a):
    z = (1./a)-1.
    v = v_from_m(m200,a)
    return prtrb.fgrowth(z,Omega_M)**0.54*5.9*v**-0.35

def v_from_m(m200,a):
    z=1./a -1.
    return (1.12*(m200/5e13)**0.3+0.53)/prtrb.fgrowth(z,Omega_M)

def Miso_r_Mdel(m,delta,a,b=0.168,scatter=False):
    #delta = del_so(delta,a)
    cc = conc(m,a)
    if(scatter):
        cc = cc* np.random.normal(1.,1./3.,len(m))
        cc = np.clip(cc,1,max(cc))
    a2 = get_A(cc)
    return 1./a2*(np.log(cc)+1./3.*np.log(delta/(3*a2*del_iso(b,a)))+1/cc*(delta/(3*del_iso(b,a)*a2))**(-1./3.)-1)

def Miso_r_Mdel_new(m_so,delta,a,b=0.168,scatter=False):
    #delta = del_so(delta,a)
    deliso = del_iso(b,a)
    cc = conc(m_so,a)
    if(scatter):
        cc = cc* np.random.normal(1.,1./3.,len(m))
        cc = np.clip(cc,1,max(cc))
    AA = get_A(cc)
    rr = ((delta/(3*deliso*AA))**(1./3.)- 2.0/(3.0*cc)) #eq 21
    MrM = 1.0/AA * (np.log(1.0+cc*rr)- (cc*rr)/(1+cc*rr)) #eq 22
    return MrM

def finite_prtcl_corr(m):
    num = m/prtcl_mass
    return (1-num**-0.6)






M = np.logspace(10,17,100)
C = conc(M,0.25)
a = get_A(C)




v=np.linspace(1,4,25)

from scipy import interpolate

#z = np.logspace(-1,1,30)-.1
a  = np.linspace(0.33,1,10)
m200 = np.logspace(10,16,10)

plt.figure()
for i in range(0,len(a)):
    plt.plot(m200,conc(m200,a[i]),label="a=%f,z=%.2f"%(a[i],dtk.z_from_a(a[i])))

plt.xscale('log')
plt.legend(loc='best')

m200_g, a_g = np.meshgrid(np.log10(m200),a)
mfof_g = Miso_r_Mdel_new(10**m200_g,200,a_g)
plt.figure()
plt.title("Mfof/M200")
for i in range(0,m200_g[:,0].size):
    plt.plot(10**m200_g[i,:],Miso_r_Mdel_new(10**m200_g[i,:],200,a_g[i,:]),'x-',label="a=%.1f, z=%.1f"%(a_g[i,0],dtk.z_from_a(a_g[i,0])))
plt.xscale('log')
plt.legend(loc='best')

m200_g = m200_g.T.flatten()
a_g    =    a_g.T.flatten()

mfof_g = np.log10(Miso_r_Mdel_new(10**m200_g,200,a_g))+m200_g
print "shape of mfof_g",mfof_g.shape

print "======"
print mfof_g
mfof_from_m200_z = interpolate.bisplrep(m200_g,a_g,mfof_g)
m200_from_mfof_z = interpolate.bisplrep(mfof_g,a_g,m200_g)


j=-1
print "a =",a_g[j]
print "true   M200: ",m200_g[j]
print "ratio: ",Miso_r_Mdel_new(10**m200_g[j],200,a_g[j])
print "true   Mfof: ",mfof_g[j],"a=",a_g[j],
print "intrpl Mfof: ",interpolate.bisplev(m200_g[j],a_g[j],mfof_from_m200_z)
#200_from_mfof_z = interpolate.griddata((mfof_g,a_g).T,
#plt.show()



mfof = np.log10(Miso_r_Mdel(m200,200,1)*finite_prtcl_corr(m200)*m200)
mfof_from_m200 = interpolate.interp1d(np.log10(m200),mfof,kind='cubic')
m200_from_mfof = interpolate.interp1d(mfof,np.log10(m200),kind='cubic')


def so_f_fof(fof_mass):
    return 10**m200_from_mfof(np.log10(fof_mass))

def fof_f_so(so_mass):
    return 10**mfof_from_m200(np.log10(so_mass))

def fof_f_so_scatter(so_mass,a,b):
    mfof = Miso_r_Mdel(so_mass,200,a,scatter=True)*so_mass
    if(b==0.2):
        mfof = mfof/finite_prtcl_corr(mfof)
    return mfof

def so_from_fof_z(fof_mass,z):
    a = 1./(z+1.0)
    log10mass = np.log10(fof_mass)
    return 10**interpolate.bisplev(log10mass,a,m200_from_mfof_z)

def fof_from_so_z(so_mass,z):
    a = 1./(z+1.0)
    log10mass = np.log10(fof_mass)
    return 10**interpolate.bisplev(log10mass,a,mfof_from_m200_z)

def rdelta_from_sod(sod_mass,delta,a=None,z=None):
    #delta here is over the over density factory,
    #typically 200
    delta = float(delta)
    rho_crit = dtk.get_rho_crit(z=z,a=a) # crit density at this redshift
    #print "rho_crit: ",rho_crit
    #m200 = (4.0*np.pi)/3.0 * r200**3 * rho_crit * delta
    rdelta = ((3.0*sod_mass)/(4.0*np.pi*rho_crit*delta))**(1./3.)
    #print "rdelta = ",rdelta
    #print "recalculat mass: %.2e"%((4.0*np.pi)/3.0*rdelta**3 * rho_crit *delta)
    return rdelta

print so_f_fof(1e13)
print so_from_fof_z(1e13,1)
print "mass interpd test: fof_real->%g, sod_real->%g sod_conv->%g fof_conv->%g"%(14,14,m200_from_mfof(14),mfof_from_m200(14))



if __name__ == "__main__":
    



# <codecell>


