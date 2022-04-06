###31-03-2022
### A script to read density and temperature arrays from LOS previously extracted and to compute 
### the corresponding values of T0 and slope (based on the idl script get_trho.txt)

import numpy as np
import os, sys
from scipy.optimize import minimize
from scipy.stats import mode


print('#***************************************************************')
print('#***************************************************************')
print('#***************************************************************')
print('#      !!!!!!!!!!!!!!!!!!! NEW RUN !!!!!!!!!!!!!!!!!!!')
print('#***************************************************************')
print('#***************************************************************')
print('#***************************************************************')

########################################## INPUT

root = "LCDM_"
labels = ['neff=-2.302']#,'neff=-2.302_Tverycold','neff=-2.302_Tveryhot']#,'s8=0.697','s8=0.967']
z_obs_list = [3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6, 5.0, 5.4]
BoxSize = 20.
num_files = 5 #number of files in which the density array is split

TEST = 'EJA' # = 'EJA' to make test on the first z-bin only
flux_path = '/scratch/rmurgia/ML_catalogues/DMNU/variousNEFF/'
#out_path = '/home/rmurgia/PF_ML/variousNEFF/'


########################generally to not to touch hereafter

nbins = 2048
nlos = 5000
nlos_per_file = nlos/num_files

#if not os.path.exists(out_path):
#  os.makedirs(out_path)

if TEST == 'EJA':

	zbins = [3.0,4.2,5.4]
	print("#You have only "+str(np.shape(zbins))+" redshift bin because it's a test!")
	print('#num of redshift bins = '+np.str(np.shape(zbins)))
	print('#num of sims = '+str(len(labels)))

else:

	zbins = z_obs_list
	print('#num of redshift bins = '+np.str(np.shape(zbins)))
	print('#num of sims = '+str(len(labels)))


############################### loop on sims and redshifts

for sim_index in range(len(labels)):   #loop on sims
	label = labels[sim_index]
	print("#*********************************************")
	print("#*START WITH model "+root+label)

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
		delta_array = np.zeros((nbins*nlos)); logdelta = np.zeros((nbins*nlos))
		
		for index in range(num_files):
			rhokerH, tempkerH1 = np.loadtxt(z_folder+root+label+"_z="+str(z_index)+"_rhotvelker"+str(index+1)+".dat", usecols=[0,2], unpack=True, comments='#')
			start = np.int(index*nbins*nlos_per_file); end = np.int((index+1)*nbins*nlos_per_file)
			density_array[start:end] = rhokerH
			temp_array[start:end] = tempkerH1

		#mean_density = np.mean(density_array)
		#print('<rhokerH> ='+str(mean_density))
		#Omega_b = 0.0457
		#h = 0.702
		#rho_crit = 1.8788*1e-26*(h**2)
		#rho_crit = 3*(h*1e2)**2/(8*np.pi)
		#rho_crit = 1.05375*1e-5*(h**2)
		#mean_density = Omega_b*rho_crit
		#print('Omega_b*rho_crit ='+str(mean_density))
	
		delta_array[:] = density_array[:]#/mean_density
		
		###OPTION 1
		#logdelta[:] = np.array(np.round(np.log10(delta_array[:]),0))
		#indices = np.argwhere(logdelta==0)
		#print(np.shape(indices))
		
		###OPTION 2
		print(np.shape(delta_array))
		logdelta = np.log10(delta_array)
		indices = np.argwhere((logdelta>-0.01) & (logdelta<0.01))
		print(np.shape(indices))
		indices = np.reshape(indices,-1)	
	
		###OPTION 3
		#vals,counts = np.unique(logdelta[indices], return_counts=True)
		#index = np.argmax(counts)
		
		T0_array = []
		for i in indices:
			T0_array.append(temp_array[i])

		T0_array = np.reshape(T0_array,-1)
		T0 = mode(T0_array)
		T0exp = mode(np.log10(T0_array))
		
		np.savetxt("./arrayT0"+str(z_index)+".txt", np.transpose([np.array(T0_array)]))
		
		###OPTION 5
		#def get_small_mode(numbers, out_mode):
		#	counts = {k:numbers.count(k) for k in set(numbers)}
		#	modes = sorted(dict(filter(lambda x: x[1] == max(counts.values()), counts.items())).keys())
		#	if out_mode=='smallest':
		#		return modes[0]
		#	elif out_mode=='largest':
		#		return modes[-1]
		#	else:
		#		return modes
	
		print("OPTION 1: T0 = "+str(T0[0]))		
		print("OPTION 2: T0 = "+str(10**T0exp[0]))
		#print("OPTION 3: T0 = "+str(T0_array[index]))
		#print("#OPTION 4: T0 = "+str(10**T0_array2[index]))
		#print("OPTION 5: T0 = "+str(get_small_mode(T0_array,'all')))

		print("#**DONE WITH z="+str(z_index))
	print("#*DONE WITH model "+root+label)
	print("#*********************************************")
