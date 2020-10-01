###01-10-2020
### A script to read tau and fluxes previously extracted to compute 
### the lyman-alpha flux power spectra

########################################## INPUT

import numpy as np
import os, sys

flux_path = '../ML_catalogues/'

out_path = '../ML_catalogues/'
if not os.path.exists(out_path):
  os.makedirs(out_path)

root = 'NCDM_'
sims = [1] #how many sims
nbins = 2048
nlos = 5000

TEST = 'EJA' # = 'EJA' to make test on 2 z-bins only

#######################################

############################### loop on sims and redshifts

if TEST == 'EJA':

	zbins = [2.2]
	print("You have only "+str(np.shape(zbins))+" redshift bin because it's a test!")
	print('num of redshift bins = '+np.str(np.shape(zbins)))
	print('num of sims = '+str(len(sims)))

else:

	zbins = []
	for i in np.linspace(2.2,6.0,num=39):
		zbins.append(i)
	print('num of redshift bins = '+np.str(np.shape(zbins)))
	print('num of sims = '+str(len(sims)))


for sim_index in sims:   #loop on sims

	sim_folder = out_path+root+str(sim_index)
	out_folder = sim_folder+'/PF'	#input/output folder
	if not os.path.exists(out_folder):
		os.makedirs(out_folder)

	for z_index in zbins:    #loop on redshifts
		
		z_index = round(z_index,2)
		z_folder = sim_folder+'/z='+str(z_index)+'/'
		
		flux_array = np.loadtxt(z_folder+root+str(sim_index)+"_z="+str(z_index)+"_tauF.dat", usecols=[1], comments='#')
		mean_flux = np.mean(flux_array)
		
		print("<F> ="+str(mean_flux))

		delta_array = np.zeros(len(flux_array))
		delta_array[:] = (flux_array[:]- mean_flux)/mean_flux 

		print(delta_array[2000:2100])



