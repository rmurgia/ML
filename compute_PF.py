###01-10-2020
### A script to read tau and fluxes previously extracted and to compute 
### the 1dimensional lyman-alpha flux power spectra

########################################## INPUT

import numpy as np
import os, sys
import matplotlib.pyplot as plt

flux_path = '../ML_catalogues/'

out_path = '../ML_catalogues/'
if not os.path.exists(out_path):
  os.makedirs(out_path)

root = 'NCDM_'
sims = [1] #how many sims
nbins = 2048
nlos = 5000
F_obs_list = [0.669181, 0.617042, 0.564612, 0.512514, 0.461362, 0.411733, 0.364155, 0.253828, 0.146033, 0.0712724]
z_obs_list = [3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6, 5.0, 5.4]

TEST = 'EJA' # = 'EJA' to make test on 2 z-bins only

############################################################################################

############################### loop on sims and redshifts

if TEST == 'EJA':

	zbins = [2.2]
	print("You have only "+str(np.shape(zbins))+" redshift bin because it's a test!")
	print('num of redshift bins = '+np.str(np.shape(zbins)))
	print('num of sims = '+str(len(sims)))

else:

	zbins = z_obs_list
	# for i in np.linspace(2.2,6.0,num=39):
	# 	zbins.append(i)
	print('num of redshift bins = '+np.str(np.shape(zbins)))
	print('num of sims = '+str(len(sims)))


for sim_index in sims:   #loop on sims

	sim_folder = out_path+root+str(sim_index)
	out_folder = sim_folder+'/PF'	#input/output folder
	if not os.path.exists(out_folder):
		os.makedirs(out_folder)

	for i,z_index in zip(range(len(z_obs_list)),zbins):    #loop on redshifts
		
		z_index = round(z_index,2)
		z_folder = sim_folder+'/z='+str(z_index)+'/'
		
		tau_array, flux_array = np.loadtxt(z_folder+root+str(sim_index)+"_z="+str(z_index)+"_tauF.dat", usecols=[0,1], unpack=True, comments='#')
		mean_flux = np.mean(flux_array)
		
		mean_flux_obs = F_obs_list[i]

		print("<F> ="+str(mean_flux))
		print("observed <F> ="+str(mean_flux_obs))

		A = np.log(mean_flux/mean_flux_obs) #compute norm. factor

		print("norm. factor A ="+str(A))

		flux_array_new = np.zeros(len(flux_array))
		delta_array = np.zeros(len(flux_array))

		flux_array_new[:] = np.exp(-A*tau_array[:])  #new flux array (normalized to the obs. flux)
		delta_array[:] = (flux_array_new[:]- mean_flux_obs)/mean_flux_obs #1D density field

		print("doing FFT...")
		##compute 1D fft of delta_array
		FFT = np.fft.fft(delta_array)
		print("..done!")

		print("compute P_F...")
		##compute 1D power spectrum
		PF = np.abs(FFT)**2
		print("..done!")

		plt.loglog(PF)
		plt.show()

