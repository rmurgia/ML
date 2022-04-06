###31-03-2022
### A script to read density and temperature arrays from LOS previously extracted and to compute 
### the corresponding values of T0 and slope (based on the idl script get_trho.txt)

import numpy as np
import os, sys
from scipy.optimize import minimize
from scipy.stats import mode


print('***************************************************************')
print('***************************************************************')
print('***************************************************************')
print('      !!!!!!!!!!!!!!!!!!! NEW RUN !!!!!!!!!!!!!!!!!!!')
print('***************************************************************')
print('***************************************************************')
print('***************************************************************')

########################################## INPUT

root = "LCDM_"
labels = ['neff=-2.302','neff=-2.302_Tverycold','neff=-2.302_Tveryhot']#,'s8=0.697','s8=0.967']
z_obs_list = [3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6, 5.0, 5.4]
BoxSize = 20.
num_files = 5 #number of files in which the density array is split

TEST = 'no' # = 'EJA' to make test on the first z-bin only
flux_path = '/scratch/rmurgia/ML_catalogues/DMNU/variousNEFF/'
#out_path = '/home/rmurgia/PF_ML/variousNEFF/'


########################generally to not to touch hereafter

nbins = 2048
nlos = 5000
nlos_per_file = nlos/num_files

#if not os.path.exists(out_path):
#  os.makedirs(out_path)

if TEST == 'EJA':

	zbins = [3.0]
	print("You have only "+str(np.shape(zbins))+" redshift bin because it's a test!")
	print('num of redshift bins = '+np.str(np.shape(zbins)))
	print('num of sims = '+str(len(labels)))

else:

	zbins = z_obs_list
	print('num of redshift bins = '+np.str(np.shape(zbins)))
	print('num of sims = '+str(len(labels)))


############################### loop on sims and redshifts

for sim_index in range(len(labels)):   #loop on sims
	label = labels[sim_index]
	print("*********************************************")
	print("*START WITH model "+root+label)

	sim_folder = flux_path+root+label
	
	#out_folder = out_path+root+label+'/T0/'	#output folder
	#if not os.path.exists(out_folder):
	#	os.makedirs(out_folder)

	for i,z_index in zip(range(len(z_obs_list)),zbins):    #loop on redshifts
		print("**START WITH z="+str(z_index))		
		
		z_index = round(z_index,2)
		z_folder = sim_folder+'/z='+str(z_index)+'/'
		
		## reading the density and temperature array
		density_array = np.zeros((nbins*nlos)); temp_array = np.zeros((nbins*nlos))
		
		for index in range(num_files):
			rhokerH, tempkerH1 = np.loadtxt(z_folder+root+label+"_z="+str(z_index)+"_rhotvelker"+str(index+1)+".dat", usecols=[0,2], unpack=True, comments='#')
			start = np.int(index*nbins*nlos_per_file); end = np.int((index+1)*nbins*nlos_per_file)
			density_array[start:end] = np.log10(rhokerH)
			temp_array[start:end] = np.log10(tempkerH1)

		nfit = 0; mx = 0; my = 0; sx = 0; sxy = 0
		
		for j in range(nbins*nlos):
			
			temp_H1 = temp_array[j]
			density = density_array[j]
			
			if density < -1.0 or density > 0.0: 
				continue
			if temp_H1 > 5.0:
				continue
			
			nfit = nfit+1
			
			mx += density
			my += temp_H1
			
			sx += density*density
			sxy += density*temp_H1
			
		beta = (sxy - mx*my/nfit)/(sx - mx*mx/nfit)
		alpha = my/nfit - beta*mx/nfit

		T0 = 10**alpha
		gamma = 1.0 + beta

		print("T0 = "+str(T0))
		print("gamma = "+str(gamma))

		print("**DONE WITH z="+str(z_index))
	print("*DONE WITH model "+root+label)
	print("*********************************************")
