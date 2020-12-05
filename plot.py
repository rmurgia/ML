import numpy as np
import matplotlib.pyplot as plt

x,y = np.loadtxt("/Users/rmurgia/Dropbox/SISSA/project/INTERPOLATION/NCDM_TEST/NCDM_TEST_1pf_z=3.0.txt", usecols=[0,1], unpack=True)

x0,y0,z0  = np.loadtxt("/Users/rmurgia/Desktop/cluster_github/ML_catalogues/NCDM_1/PF/PF_z3.dat", usecols=[0,1,2], unpack=True)

# plt.loglog(x,y, label="old")
# plt.loglog(x0,y0, label="new, dim")
# plt.loglog(x0,z0, label="new, adim")
plt.xscale('log')
plt.plot(x0,(y-y0)/y)
plt.ylim(-0.01,0.01)
plt.legend()
plt.show()
