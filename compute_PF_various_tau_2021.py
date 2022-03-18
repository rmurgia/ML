###01-10-2020
### A script to read tau and fluxes previously extracted and to compute 
### the 1dimensional lyman-alpha flux power spectra

import numpy as np
import os, sys
from scipy.optimize import minimize

########################################## INPUT

#root = 'lcd'
#labels = ['m']

#root = 'PBHs_'
#labels = ['1e1-5_NEW','1e2-2_NEW','1e2-3_NEW', '1e2', '1e2-5','1e2-7','1e3-5','1e1','1e2-4','1e2-6','1e3','1e4']
#labels = ['1e3']
root = "LCDM_neff="
labels = ['-2.302','-2.453','-2.583']

F_obs_list = [0.669181, 0.617042, 0.564612, 0.512514, 0.461362, 0.411733, 0.364155, 0.253828, 0.146033, 0.0712724]
z_obs_list = [3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6, 5.0, 5.4]
#F_obs_list = [0.364155, 0.253828, 0.146033, 0.0712724]
# z_obs_list = [4.2, 4.6, 5.0, 5.4]

BoxSize = 20.

## to rescale wrt to mean fluxes different from the observed ones
NONST_RESCALE = 'YES'
print("NONST_RESCALE = "+NONST_RESCALE)

#FACTORS = [0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4]
FACTORS = [0.6, 0.8, 1.2, 1.4]
#FACTORS = [1.4]

for FACTOR in FACTORS:

	if NONST_RESCALE == 'YES':
	    F_obs_list = np.array(F_obs_list)
	    tau_obs_list = -np.log(F_obs_list)
	    tau_obs_list_rescaled = FACTOR*tau_obs_list
	    F_obs_list = np.exp(-tau_obs_list_rescaled)
	
	print("********")
	
	TEST = 'no' # = 'EJA' to make test on the first z-bin only
	flux_path = '/scratch/rmurgia/ML_catalogues/DMNU/variousNEFF/'
	out_path = '/home/rmurgia/PF_ML/variousNEFF/'
	#flux_path = '../ML_catalogues/'
	#out_path = '../ML_catalogues/'
	if not os.path.exists(out_path):
	  os.makedirs(out_path)
	
	
	########################generally to not to touch hereafter
	nbins = 2048
	nlos = 5000
	PF_matrix = np.zeros((nlos,nbins))
	PF = np.zeros(nbins)
	kF = 2.0*np.pi/BoxSize #fundamental frequency  
	middle = nbins/2  
	kN = middle*kF #Nyquist frequency
	spatial_step = 0.5/kN
	
	def func_A(A, tau_sim_list, meanF_obs, meanF_sim):
		factor = meanF_sim/meanF_obs
		numerator = np.sum(np.exp(-tau_sim_list))
		denominator = np.sum(np.exp(-tau_sim_list*A))
		result = np.abs(factor - (numerator/denominator))
		return result
	###################################################################
	
	
	if TEST == 'EJA':
	
		zbins = [3.0]
		print("You have only "+str(np.shape(zbins))+" redshift bin because it's a test!")
		print('num of redshift bins = '+np.str(np.shape(zbins)))
		print('num of sims = '+str(len(labels)))
	
	############################### loop on sims and redshifts
	else:
	
		zbins = z_obs_list
		print('num of redshift bins = '+np.str(np.shape(zbins)))
		print('num of sims = '+str(len(labels)))
	
	
	for sim_index in range(len(labels)):   #loop on sims
	
		label = labels[sim_index]
		sim_folder = flux_path+root+label
		out_folder = out_path+root+label+'/PF/'	#input/output folder
		if not os.path.exists(out_folder):
			os.makedirs(out_folder)
	
		for i,z_index in zip(range(len(z_obs_list)),zbins):    #loop on redshifts
			
			z_index = round(z_index,2)
			z_folder = sim_folder+'/z='+str(z_index)+'/'
			
			tau_array, flux_array = np.loadtxt(z_folder+root+label+"_z="+str(z_index)+"_tauF.dat", usecols=[0,1], unpack=True, comments='#')
			mean_flux = np.mean(flux_array)
			
			mean_flux_obs = F_obs_list[i]
	
			print("<F> ="+str(mean_flux))
			print("observed <F> ="+str(mean_flux_obs))
	
			if mean_flux < mean_flux_obs:
				A_list = np.linspace(0.01,1., num=100)
				print("A < 1")
			
			elif mean_flux > mean_flux_obs:
				A_list = np.linspace(1.,10.0, num=100)
				print("A > 1")
			
			y = []
			for i,A in zip(range(len(A_list)),A_list):
				y.append(func_A(A, tau_array, mean_flux_obs, mean_flux))
			
			print("minimum difference before optimization is "+str(min(y)))
			
			res = minimize(func_A, min(y), args=(tau_array, mean_flux_obs, mean_flux), tol=1e-5)
			best_A = res.x[0]
			print("Start with optimization...")
			print(res)
			print("..Optimization done!")
			
			print("norm. factor A ="+str(best_A))
	
			flux_array_new = np.zeros(len(flux_array))
			delta_array = np.zeros(len(flux_array))
	
			flux_array_new[:] = np.exp(-tau_array[:]*best_A)  #new flux array (normalized to the obs. flux)
			
			mean_flux_new = np.mean(flux_array_new)
			print("<F>_NEW ="+str(mean_flux_new))
			delta_array[:] = (flux_array_new[:] - mean_flux_new)/mean_flux_new #1D density field
			
			print("len(deltas) = "+str(len(delta_array)))
			
			#compute FFT and P_flux per each los
			for i in range(nlos):		
				deltas = delta_array[nbins*i:nbins*(i+1)]
				
				if len(deltas) != nbins:
					print("ERROR!!"); break
	
				##compute 1D fft of delta_array
				FFT = np.fft.fft(deltas)/len(deltas)
				
				##compute 1D power spectrum
				PF_matrix[i][:] = np.abs(FFT)**2
				
			PF = PF_matrix[0][:]
			freqs = np.fft.fftfreq(deltas.size, spatial_step)
			idx = np.argsort(freqs)
	
			print(np.shape(freqs),np.shape(PF))
			
			## average over the nlos different P_flux
			to_be_avg = np.zeros(nlos)
			for i in range(nbins):
			 	to_be_avg = PF_matrix[:,i]
			 	PF[i] = np.mean(to_be_avg)
	
			end = int(0.5*nbins+1)
	
			##converting from dimension-less power spectrum to P_F(k) in (h/Mpc) (WE ARE IN 1-DIM)
			freqs_final = np.abs(freqs[1:end])
			PF_final = np.zeros(len(freqs_final))
			#PF_final[:] = 2*np.pi*PF[1:end]/freqs_final[:]
			PF_final[:] = PF[1:end]
			
			if NONST_RESCALE == 'YES':
				if label == '1e1-5_NEW':
					label = '1e1-5'
				elif label == '1e2-2_NEW':
					label = '1e2-2'
				elif label == '1e2-3_NEW':
					label = '1e2-3'
				np.savetxt(out_folder+"PF_"+root+label+"_F"+str(FACTOR)+"_z"+str(z_index)+".dat",np.transpose([freqs_final,PF_final]))
				if label == '1e1-5':
					label = '1e1-5_NEW'
				elif label == '1e2-2':
					label = '1e2-2_NEW'
				elif label == '1e2-3':
					label = '1e2-3_NEW'
			else:
				np.savetxt(out_folder+"PF_"+root+label+"_z"+str(z_index)+".dat",np.transpose([freqs_final,PF_final]))
			print("**DONE WITH z="+str(z_index))
		print("*DONE WITH model "+root+label)
	
