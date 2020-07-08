#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global_vars.h"
#include "proto.h"


double ztime[1],omegam[1],omegal[1],omegab[1],h100[1],box100[1],Xh[1];
int nbins[1],nlos[1];
  
double *rhoker_H,*posaxis,*velaxis;
double *rhoker_H1,*velker_H1,*tempker_H1,*tau_H1;
double *xlos,*ylos,*zlos;
int *ilos;


double ztime_file;
char *path;

/********************************************************************************/

int main(int argc, char **argv)
{

  int i;

  printf("\n");
  
  if(argc != 2)
    {
      fprintf(stderr, "\n\nIncorrect argument(s).  Specify:\n\n");
      fprintf(stderr, "<path>     (path to output files)\n");
      fprintf(stderr, "Usage: ./LosExtract <path>\n");
      fprintf(stderr, "\n\n");
      exit(0);
    }
  
  path = argv[1];
  
   ztime_file = 2.200;
   //ztime_file = 2.000;
  //ztime_file = 5.200;

   for(i=0; i<40; i++)
  //for(i=0; i<101; i++)
  //  for(i=0; i<69; i++)
    {
      printf("\nReading los output at z=%.3f:\n",ztime_file);
      read_los();
      
      
      printf("Computing optical depths...\n");
      compute_absorption();
      printf("Done.\n\n");
      
      write_tau();
      
      free_los_memory();

      ztime_file += 0.1;
    }
  
  return(0);
}

/********************************************************************************/

void compute_absorption()
{
  
  double T0,T1,T2;
  double atime,vmax,vmax2,rscale,escale,drbin,Hz;
  double H0,rhoc,critH;
  double pcdone,pccount=0.0;
  double sigma_Lya_H1,k1_H1,k2_H1,k3_H1;
  double profile_H1,vdiff_H1;
  double u_H1[nbins[0]],binv_H1[nbins[0]],b2inv_H1[nbins[0]];
  double aa_H1[nbins[0]],tau_delta_H1[nbins[0]];
  
  int i,j,iproc,convol_index,pixel_index;
  
  atime  = 1.0/(1.0+ztime[0]);
  rscale = (KPC*atime)/h100[0];        /* comoving kpc/h to cm */
  escale = 1.0e10;                     /* (km s^-1)^2 to (cm s^-1)^2 */
  drbin  = box100[0]/(double)nbins[0]; /* comoving kpc/h */
  Hz     = 100.0*h100[0]*sqrt(omegam[0]/(atime*atime*atime)+omegal[0]); /* km s^-1 Mpc^-1 */
  vmax   = box100[0]*Hz*rscale/MPC; /* box size km s^-1 */
  vmax2  = 0.5*vmax;
  H0     = 1.0e7/MPC; /* 100 km s^-1 Mpc^-1 in cgs */ 
  rhoc   = 3.0*(H0*h100[0])*(H0*h100[0])/(8.0*PI*GRAVITY); /* g cm^-3 */
  critH  = rhoc*omegab[0]*Xh[0]/(atime*atime*atime); /* g cm^-3*/
  
  
  /* HI Lyman-alpha */
  sigma_Lya_H1 = sqrt(3.0*PI*SIGMA_T/8.0)*LAMBDA_LYA_H1*FOSC_LYA; /* cm^2 */
  k1_H1        = 2.0*BOLTZMANN/(HMASS*AMU);
  k2_H1        = GAMMA_LYA_H1*LAMBDA_LYA_H1/(4.0*PI);
  k3_H1        = sigma_Lya_H1*C*rscale*drbin*critH/(sqrt(PI)*HMASS*AMU);

  for(iproc=0; iproc<nlos[0]; iproc++)
    {
      for(j=0; j<nbins[0]; j++)
  	{
  	  convol_index =  j + nbins[0]*iproc;
	  
  	  /* HI Lyman-alpha */
  	  u_H1[j]         = velaxis[j] + velker_H1[convol_index]; /* km s^-1 */
  	  binv_H1[j]      = 1.0/sqrt(k1_H1*tempker_H1[convol_index]); /* cm s^-1 */
  	  b2inv_H1[j]     = escale*binv_H1[j]*binv_H1[j];  /* (km s^-1)^-2 */
  	  aa_H1[j]        = k2_H1*binv_H1[j];
  	  tau_delta_H1[j] = k3_H1*rhoker_H[convol_index]*rhoker_H1[convol_index]*binv_H1[j];
  	}
      
      /* Voigt profile: Tepper-Garcia, 2006, MNRAS, 369, 2025 */
      for(i=0; i<nbins[0]; i++)
  	{
  	  pixel_index =  i + nbins[0]*iproc;
	  
  	  for(j=0; j<nbins[0]; j++)
  	    {
	      /* HI Lyman-alpha */
  	      vdiff_H1 = fabs(velaxis[i] - u_H1[j]);
  	      if (vdiff_H1 > vmax2) vdiff_H1 = vmax - vdiff_H1;
	      
  	      T0 = vdiff_H1*vdiff_H1*b2inv_H1[j];
  	      T1 = exp(-T0);
  	      T2 = 1.5/T0;
	      
  	      profile_H1 = (T0 < 1.0e-6)
  		? T1 : T1-aa_H1[j]/sqrt(PI)/T0*(T1*T1*(4.0*T0*T0+7.0*T0+4.0+T2)-T2-1.0);
	      
  	      tau_H1[pixel_index] += tau_delta_H1[j]*profile_H1;
	    }
  	}
      
      pcdone = 100.0*(float)iproc/((float)nlos[0]-1);
      if(pcdone >= pccount)
  	{
  	  printf("%3.2f%%\n",pcdone);
  	  pccount += 10.0;
  	}
    } 
}
  
/********************************************************************************/

void read_los()
{
  char fname[400];
  FILE *input;

  sprintf(fname, "%s/los2048_n5000_z%.3f.dat",path,ztime_file);
  if(!(input=fopen(fname,"rb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(0);
    }
  
  fread(ztime,sizeof(double),1,input);
  fread(omegam,sizeof(double),1,input);
  fread(omegal,sizeof(double),1,input);
  fread(omegab,sizeof(double),1,input);
  fread(h100,sizeof(double),1,input);
  fread(box100,sizeof(double),1,input);
  fread(Xh,sizeof(double),1,input); /* hydrogen fraction by mass */
  fread(nbins,sizeof(int),1,input);
  fread(nlos,sizeof(int),1,input);
  
  allocate_los_memory();
  
  fread(ilos,sizeof(int),nlos[0],input); /* LOS axis */
  fread(xlos,sizeof(double),nlos[0],input); /* LOS positions, comoving kpc/h */
  fread(ylos,sizeof(double),nlos[0],input);
  fread(zlos,sizeof(double),nlos[0],input);
  fread(posaxis,sizeof(double),nbins[0],input); /* pixel positions, comoving kpc/h */
  fread(velaxis,sizeof(double),nbins[0],input); /* pixel positions, km s^-1 */
  fread(rhoker_H,sizeof(double),nbins[0]*nlos[0],input); /* gas overdensity, Delta=rho/rho_crit */
  
  fread(rhoker_H1,sizeof(double),nbins[0]*nlos[0],input); /* n_HI/n_H */
  fread(tempker_H1,sizeof(double),nbins[0]*nlos[0],input); /* T [K], HI weighted */
  fread(velker_H1,sizeof(double),nbins[0]*nlos[0],input); /* v_pec [km s^-1], HI weighted */

  fclose(input);
}

/********************************************************************************/

void write_tau()
{
  char fname[400];
  FILE *output;

  sprintf(fname, "%s/tau2048_n5000_z%.3f.dat",path,ztime_file);
  if(!(output=fopen(fname,"wb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(0);
    }
  fwrite(tau_H1,sizeof(double),nbins[0]*nlos[0],output);
  fclose(output);
}

/********************************************************************************/

void allocate_los_memory()
{  
  ilos      = (int *)calloc(nlos[0], sizeof(int));
  if(NULL==ilos){free(ilos);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}
  xlos      = (double *)calloc(nlos[0], sizeof(double));
  if(NULL==xlos){free(xlos);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}
  ylos      = (double *)calloc(nlos[0], sizeof(double));
  if(NULL==ylos){free(ylos);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}
  zlos      = (double *)calloc(nlos[0], sizeof(double));
  if(NULL==zlos){free(zlos);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}
  posaxis      = (double *)calloc(nbins[0], sizeof(double));
  if(NULL==posaxis){free(posaxis);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}
  velaxis      = (double *)calloc(nbins[0], sizeof(double));
  if(NULL==posaxis){free(posaxis);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}
  rhoker_H     = (double *)calloc(nlos[0]*nbins[0], sizeof(double));
  if(NULL==rhoker_H){free(rhoker_H);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}
  
  
  rhoker_H1    = (double *)calloc(nlos[0]*nbins[0], sizeof(double));
  if(NULL==rhoker_H1){free(rhoker_H1);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}
  velker_H1    = (double *)calloc(nlos[0]*nbins[0], sizeof(double));
  if(NULL==velker_H1){free(velker_H1);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}
  tempker_H1   = (double *)calloc(nlos[0]*nbins[0], sizeof(double));
  if(NULL==tempker_H1){free(tempker_H1);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}
  tau_H1       = (double *)calloc(nlos[0]*nbins[0], sizeof(double));
  if(NULL==tau_H1){free(tau_H1);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}
}

/********************************************************************************/

void free_los_memory()
{  
  free(ilos);
  free(xlos);
  free(ylos);
  free(zlos);
  free(posaxis);
  free(velaxis);
  free(rhoker_H);
  
  free(rhoker_H1);
  free(velker_H1);
  free(tempker_H1);
  free(tau_H1);
}

/********************************************************************************/
