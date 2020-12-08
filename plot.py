import numpy as np
import matplotlib.pyplot as plt

# [3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.6,5.0]
for i in [4.2]:#4.2,4.6,5.0,5.4]:
	
	x,y = np.loadtxt("/Users/rmurgia/Dropbox/SISSA/project/INTERPOLATION/NCDM_TEST/NCDM_TEST_1pf_z="+str(i)+".txt", usecols=[0,1], unpack=True)

	# x0,y0  = np.loadtxt("/Users/rmurgia/Desktop/cluster_github/ML_catalogues/NCDM_1/PF/PF_NCDM__1_z"+str(i)+".dat", usecols=[0,1], unpack=True)
	x0,z0  = np.loadtxt("/Users/rmurgia/Desktop/cluster_github/ML_catalogues/NCDM_1/PF/PF_NCDM_1_z"+str(i)+".dat", usecols=[0,1],unpack=True)
	

	plt.loglog(x,y, label="old")
	# plt.loglog(x0,y0, label="new")
	plt.loglog(x0,z0, label="new new")
	# y = (2*np.pi)*y/x
	plt.title("z = "+str(i))
	# plt.xscale('log')
	# plt.xlim(min(x),1e2)
	# plt.plot(x,(y-z0)/y)
	# plt.ylim(-0.2,0.2)
	
	plt.legend()
	plt.show()
