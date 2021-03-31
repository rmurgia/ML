from ctypes import *
import numpy as np

import matplotlib.pyplot as plt

import struct
import os, sys

###################################### DEFINITIONS

def read_chunks(f, length):
  while True:
    data = f.read(length)
    if not data: break
    yield data

class Read_Tau(Structure):
  _fields_ = [('tauH1', c_double)]


def len_all(ilos,ylos,zlos,posaxis,velaxis,rhokerH,rhokerH1,tempkerH1,velkerH1): 
	return(len(ilos),len(xlos),len(ylos),len(zlos),len(posaxis),len(velaxis),np.shape(rhokerH),np.shape(rhokerH1),np.shape(tempkerH1),np.shape(velkerH1))
########################################## INPUT

sim_path = '/scratch/gscelfo/PBHs_sims_z199_box20/'

out_path = '/scratch/rmurgia/ML_catalogues/'
if not os.path.exists(out_path):
  os.makedirs(out_path)

root = 'lcd'
labels = ['m']
sim_num = len(labels) #how many sims

steps = 5   #number of files where you want to store "ker" quantities
TEST = 'EJA' # = 'EJA' to make test on 2 z-bins only

############################### loop on sims and redshifts

if TEST == 'EJA':

	zbins = [2.2, 6.0]
	print("You have only "+str(np.shape(zbins))+" redshift bins because it's a test!")
	print('num of reshift bins = '+np.str(np.shape(zbins)))
	print('num of sims = '+str(sim_num))

else:

	zbins = []
  	for i in np.linspace(2.1,6.0,num=40):
    		zbins.append(i)
 	zbins.append(9.0)
  	print('num of reshift bins = '+np.str(np.shape(zbins)))
  	print('num of sims = '+str(sim_num))


for sim_index in range(sim_num):   #loop on sims

	label = labels[sim_index]	
	los_path = sim_path+root+label+'/los/'
	print("I'm starting with the "+label+" sim")

	sim_out_folder = out_path+root+str(label)+'/'	#making output folder folder for the sim
	if not os.path.exists(sim_out_folder):
    		os.makedirs(sim_out_folder)

	for z_index in zbins:    #loop on redshifts
		
		z_index = round(z_index,2)
		z_out_folder = sim_out_folder+'z='+str(z_index)+'/' 	#making output folder folder for the z-bin
		if not os.path.exists(z_out_folder):
      			os.makedirs(z_out_folder)

    ####################### stuff needed to read LOS

	  	los_file = los_path+'los2048_n5000_z'+str(z_index)+'00.dat'		

		struct_fmt = '=dddddddii5000i5000d5000d5000d2048d2048d10240000d10240000d10240000d10240000d'
		struct_len = struct.calcsize(struct_fmt)
		struct_unpack = struct.Struct(struct_fmt).unpack_from

		with open(los_file, "rb") as f:   #READING FILE
    			results = [struct_unpack(chunk) for chunk in read_chunks(f, struct_len)]

		print("los has been read (& blue):")
		print("LOS file elements = "+str(np.size(results)))
		print( np.shape(results))
	
		#################### read header and compute offset

		header = 9
		nbins = int(results[0][7])
		nlos = int(results[0][8])

		h_labels=['ztime','om_out','ol_out','ob_out','h_out','box_out','xh_out']

		print("Header Read (& Blue)")

		print("n bins = "+str(nbins))
		print("n los = "+str(nlos))

		for i in range(header-2):
			print(h_labels[i]+" = "+str(results[0][i])

		nbins_range = np.linspace(1,nbins, num=nbins)
		nlos_range = np.linspace(1,nlos, num=nlos)


		if z_index == 2.0 or z_index == 2.1 or z_index == 9.0:   ##ONLY LOS


      			### READ (AND PLOT) LOS
      			ilos = []; xlos = []; ylos = []; zlos = []; xlos_arr = []; ylos_arr = []; zlos_arr = []
      			posaxis = []; velaxis = []
      			rhokerH = []; rhokerH1 = []; tempkerH1 = []; velkerH1 = []
      			for i in range(4*nlos+6):
        
        			if i < 4: #ilos, xlos, ylos, zlos (+offset due to header)
          				start = header
          				y = results[0][start+nlos*i:start+nlos*(i+1)]

          				if i == 0:
            					ilos.append(y)
            					ilos = np.reshape(ilos,-1)
          				elif i == 1:
            					xlos.append(y)
            					xlos = np.reshape(xlos,-1)
          				elif i == 2:
            					ylos.append(y)
            					ylos = np.reshape(ylos,-1)
          				elif i == 3:
            					zlos.append(y)
            					zlos = np.reshape(zlos,-1)

        			elif 4 <= i < 6: #posaxis, velaxis
          				start1 = header + 4*nlos
          				y =  results[0][start1+nbins*(i-4):start1+nbins*(i-3)]

         				if i == 4:
            					posaxis.append(y)
            					posaxis = np.reshape(posaxis,-1)
          				elif i == 5:
            					velaxis.append(y)
            					velaxis = np.reshape(velaxis,-1)

       		 		elif 6 <= i < nlos+6: #rhokerH
          				start2 = header + 4*nlos + 2*nbins
          				y =  results[0][start2+nbins*(i-6):start2+nbins*(i-5)]
          				rhokerH = np.concatenate((rhokerH,y))
          
        			elif nlos+6 <= i < (2*nlos)+6: #rhokerH1
          				offset = header + 4*nlos + 2*nbins
          				start3 = offset + nbins*nlos
          				y =  results[0][start3+nbins*(i-6-nlos):start3+nbins*(i-5-nlos)]
          				rhokerH1 = np.concatenate((rhokerH1,y))
          
        			elif 2*nlos+6 <= i < (3*nlos)+6: #tempkerH1
          				start3 = offset + 2*nbins*nlos
          				y =  results[0][start3+nbins*(i-6-2*nlos):start3+nbins*(i-5-2*nlos)]
          				tempkerH1 = np.concatenate((tempkerH1,y))
          
        			elif 3*nlos+6 <= i < (4*nlos)+6: #velkerH1
          				start3 = offset + 3*nbins*nlos
          				y =  results[0][start3+nbins*(i-6-3*nlos):start3+nbins*(i-5-3*nlos)]
          				velkerH1 = np.concatenate((velkerH1,y))
          
            		## save results to file (one per each redshift bin)
			fout1 = z_out_folder+root+label+'_z='+str(z_index)+'_xyzlos.dat'
      			fout2 = z_out_folder+root+label+'_z='+str(z_index)+'_pvaxis.dat'
      			fout3nodat = z_out_folder+root+label+'_z='+str(z_index)+'_rhotvelker'
      			fout3 = z_out_folder+root+label+'_z='+str(z_index)+'_rhotvelker.dat'

      			FH = "n bins = "+str(nbins)+" \n n los = "+str(nlos)+" \n "+str(h_labels[:7])+" \n "+str(results[0][:7])+" \n"
      			h1 = "ilos, xlos, ylos, zlos"
      			h2 = "posaxis, velaxis"
      			h3 = "rhoker_H, rhoker_H1 tempker_H1, velker_H1"

      			for i in range(len(ilos)):
	     			if ilos[i] == 1:
				  xlos_arr.append(0.)
				  ylos_arr.append(ylos[i])
				  zlos_arr.append(zlos[i])
        			elif ilos[i] == 2:
				  xlos_arr.append(xlos[i])
				  ylos_arr.append(0.)
				  zlos_arr.append(zlos[i])
        			elif ilos[i] == 3:
				  xlos_arr.append(xlos[i])
				  ylos_arr.append(ylos[i])
				  zlos_arr.append(0.)

      			for j in range(1,4):
        			if j == 1:
          				np.savetxt(fout1, np.transpose([xlos_arr,ylos_arr,zlos_arr]), header = FH+h1)
          				print("> File %s saved: " %fout1)
        			elif j == 2:
          				np.savetxt(fout2, np.transpose([posaxis,velaxis]), header = FH+h2)
          				print("> File %s saved: " %fout2)
          				print('steps = '+str(steps))
        			elif j == 3:
          				l_step = (nlos/steps)*nbins
          				for f_index in range(0,steps):
            					np.savetxt(fout3nodat+str(f_index+1)+'.dat', np.transpose([rhokerH[f_index*l_step:(f_index+1)*l_step], rhokerH1[f_index*l_step:(f_index+1)*l_step], tempkerH1[f_index*l_step:(f_index+1)*l_step], velkerH1[f_index*l_step:(f_index+1)*l_step]]), header = FH+h3)
            				print("> File %s saved: " %fout3)
           
		else:      ##BOTH LOS AND TAU

    			################ stuff needed to read TAU

      			tau_file = los_path+'tau2048_n5000_z'+str(z_index)+'00.dat'
           
			tau_arr = []
      			flux_arr = []
 
      			with open(tau_file, 'rb') as file:  #READING FILE
        			x = Read_Tau()
        			while file.readinto(x) == sizeof(x):
          				tau_arr.append(x.tauH1)
          				flux_arr.append(np.exp(-x.tauH1))

      			print("tau has been read (& blue):")
      			print("TAU file elements = "+str(np.size(tau_arr)))
      			print(np.shape(tau_arr))


      			### READ (AND PLOT) LOS
      			ilos = []; xlos = []; ylos = []; zlos = []; xlos_arr = []; ylos_arr = []; zlos_arr = []
      			posaxis = []; velaxis = []
      			rhokerH = []; rhokerH1 = []; tempkerH1 = []; velkerH1 = []
      			for i in range(4*nlos+6):
        
        			if i < 4: #ilos, xlos, ylos, zlos (+offset due to header)
          				start = header
          				y = results[0][start+nlos*i:start+nlos*(i+1)]

          				if i == 0:
            					ilos.append(y)
            					ilos = np.reshape(ilos,-1)
          				elif i == 1:
            					xlos.append(y)
            					xlos = np.reshape(xlos,-1)
          				elif i == 2:
            					ylos.append(y)
            					ylos = np.reshape(ylos,-1)
          				elif i == 3:
            					zlos.append(y)
            					zlos = np.reshape(zlos,-1)

        			elif 4 <= i < 6: #posaxis, velaxis
          				start1 = header + 4*nlos
          				y =  results[0][start1+nbins*(i-4):start1+nbins*(i-3)]

          				if i == 4:
            					posaxis.append(y)
           					posaxis = np.reshape(posaxis,-1)
          				elif i == 5:
           	 				velaxis.append(y)
            					velaxis = np.reshape(velaxis,-1)

        			elif 6 <= i < nlos+6: #rhokerH
          				start2 = header + 4*nlos + 2*nbins
          				y =  results[0][start2+nbins*(i-6):start2+nbins*(i-5)]
          				rhokerH = np.concatenate((rhokerH,y))
          
        			elif nlos+6 <= i < (2*nlos)+6: #rhokerH1
          				offset = header + 4*nlos + 2*nbins
          				start3 = offset + nbins*nlos
        		  		y =  results[0][start3+nbins*(i-6-nlos):start3+nbins*(i-5-nlos)]
          				rhokerH1 = np.concatenate((rhokerH1,y))
          
        			elif 2*nlos+6 <= i < (3*nlos)+6: #tempkerH1
          				start3 = offset + 2*nbins*nlos
          				y =  results[0][start3+nbins*(i-6-2*nlos):start3+nbins*(i-5-2*nlos)]
          				tempkerH1 = np.concatenate((tempkerH1,y))
          
        			elif 3*nlos+6 <= i < (4*nlos)+6: #velkerH1
          				start3 = offset + 3*nbins*nlos
          				y =  results[0][start3+nbins*(i-6-3*nlos):start3+nbins*(i-5-3*nlos)]
          				velkerH1 = np.concatenate((velkerH1,y))
          
            		## save results to file (one per each redshift bin)
      			fout1 = z_out_folder+root+label+'_z='+str(z_index)+'_xyzlos.dat'
      			fout2 = z_out_folder+root+label+'_z='+str(z_index)+'_pvaxis.dat'
			fout3nodat = z_out_folder+root+label+'_z='+str(z_index)+'_rhotvelker'
      			fout3 = z_out_folder+root+label+'_z='+str(z_index)+'_rhotvelker.dat'

      			FH = "n bins = "+str(nbins)+" \n n los = "+str(nlos)+" \n "+str(h_labels[:7])+" \n "+str(results[0][:7])+" \n"
      			h1 = "ilos, xlos, ylos, zlos"
      			h2 = "posaxis, velaxis"
      			h3 = "rhoker_H, rhoker_H1 tempker_H1, velker_H1"

      			for i in range(len(ilos)):
        			if ilos[i] == 1:
          				xlos_arr.append(0.)
          				ylos_arr.append(ylos[i])
          				zlos_arr.append(zlos[i])
        			elif ilos[i] == 2:
          				xlos_arr.append(xlos[i])
          				ylos_arr.append(0.)
          				zlos_arr.append(zlos[i])
        			elif ilos[i] == 3:
          				xlos_arr.append(xlos[i])
          				ylos_arr.append(ylos[i])
          				zlos_arr.append(0.)

      			for j in range(1,4):
        			if j == 1:
          				np.savetxt(fout1, np.transpose([xlos_arr,ylos_arr,zlos_arr]), header = FH+h1)
          				print("> File %s saved: " %fout1)
        			elif j == 2:
          				np.savetxt(fout2, np.transpose([posaxis,velaxis]), header = FH+h2)
          				print("> File %s saved: " %fout2)
          				print('steps = '+str(steps))
        			elif j == 3:
          				l_step = (nlos/steps)*nbins
          				for f_index in range(0,steps):
            					np.savetxt(fout3nodat+str(f_index+1)+'.dat', np.transpose([rhokerH[f_index*l_step:(f_index+1)*l_step], rhokerH1[f_index*l_step:(f_index+1)*l_step], tempkerH1[f_index*l_step:(f_index+1)*l_step], velkerH1[f_index*l_step:(f_index+1)*l_step]]), header = FH+h3)
            			print("> File %s saved: " %fout3)			

    			### READ (AND PLOT) TAU
			fout = z_out_folder+root+label+'_z='+str(z_index)+'_tauF.dat'
    			FH = "n bins = "+str(nbins)+" \n n los = "+str(nlos)+" \n "+str(h_labels[:7])+" \n "+str(results[0][:7])+" \n tau_arr, flux_arr"
			np.savetxt(fout, np.transpose([tau_arr,flux_arr]), header = FH)
      			print("> File %s saved: " %fout)

		print("z = "+str(z_index)+" for sim "+root+label+" read (& blue)!")

	print("sim "+root+label+" completely read (& blue)!")
