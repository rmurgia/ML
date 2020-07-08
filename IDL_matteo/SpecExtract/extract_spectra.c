#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global_vars.h"
#include "parameters.h"
#include "proto.h"



double *rhoker_H,*posaxis,*velaxis;
double *rhoker_H1,*velker_H1,*tempker_H1,*tau_H1;
double *xlos,*ylos,*zlos;
int *ilos;

double vmax,vmax2,drbin;


/********************************************************************************/

/* Perform the SPH interpolation along the sight-lines */
void SPH_interpolation()
{
  double Hz,mu,nelec;
  double drinv,box2;
  double dscale,rscale,vscale,mscale,escale;
  double xx,yy,zz,hh,h2,h4,dr,dr2;
  double hinv2,hinv3,drmax,rgrid;
  double dist2,q,kernel,velker,tempker;
  double pcdone,pccount=0.0;
  
  int i,iproc,iz=0,ioff,j,iiz;
  int pixel_index;
  
  
 
  Hz = 100.0*h100*sqrt(omega0/(atime*atime*atime)+omegaL); /* km s^-1 Mpc^-1 */
  
  
  /* Conversion factors from internal units */
  rscale = (KPC*atime)/h100;              /* convert length to cm */
  vscale = sqrt(atime);                   /* convert velocity to kms^-1 */
  mscale = (1.0e10*SOLAR_MASS)/h100;      /* convert mass to g */
  escale = 1.0e10;                        /* convert energy/unit mass from(km/s)^2 to (cm/s)^2 */
  dscale = mscale/(rscale*rscale*rscale); /* convert density to g cm^-3 */
  
  /* Box scales in internal units */
  drbin = box100/(double)NBINS; /* bin size */
  drinv = 1.0/drbin;
  box2  = 0.5*box100;
  
  /* Box scales in km s^-1 */
  vmax  = box100*Hz*rscale/MPC; /* box size */
  vmax2 = 0.5*vmax;
  
  
  for(i=0;i<Ntype;i++)
    {
      nelec = P[i].Ne;
      mu      = 1.0/(XH*(0.75+nelec) + 0.25);  /* Mean molecular weight */
      P[i].U *= ((GAMMA-1.0)*mu*HMASS*AMU)/BOLTZMANN; /* K/escale */
      P[i].h *= 0.5; /* Note the factor of 0.5 for this kernel definition */
    }
  
  
  /* SPH interpolation onto sight-lines */
  for(iproc=0; iproc<NUMLOS; iproc++) 
    {
      for(i=0; i<Ntype; i++)
	{
	  xx = P[i].Pos[0];
	  yy = P[i].Pos[1];
	  zz = P[i].Pos[2];
	  
	  hh = P[i].h;
	  
	  /* distance to projection axis */

	  dr = (ilos[iproc] == 1) ? fabs(yy-ylos[iproc]) : fabs(xx-xlos[iproc]);
	  
	  if (dr > box2) dr = box100 - dr; 
	  
	  /* particle is within 2h of the axis */
	  if (dr <= 2.0*hh) 
	    {
	      dr2 = dr*dr;
	      
	      /* distance to projection axis */	
	      dr = (ilos[iproc] == 3) ? fabs(yy-ylos[iproc]) : fabs(zz-zlos[iproc]); 
	      
	      if (dr > box2) dr = box100 - dr; 
              
	      dr2 += dr*dr; /* dx^2+dy^2 */
	      
	      h2 = hh*hh; 
	      h4 = 4.0*h2; 
	      
	      /* particle is within 2h of the axis */
	      if (dr2 <= h4)
		{
		  hinv2 = 1.0   / h2; /* 1/h^2 */
		  hinv3 = hinv2 / hh; /* 1/h^3 */
		  
		  /* Central pixel to contribute to */
		  if (ilos[iproc] == 1)
		    iz = xx * drinv + 0.5;
		  else if (ilos[iproc] == 2) 
		    iz = yy * drinv + 0.5;
		  else if (ilos[iproc] == 3)
		    iz = zz * drinv + 0.5;
		  
		  
		  /* max distance along line-of-sight axis from central
		     pixel lying within 2h of contributing particle */
		  drmax = sqrt(fabs(h4 - dr2)); 
		  ioff  = (int)(drmax * drinv);
		  
		  /* Loop over contributing pixels */
		  for(iiz=iz-ioff; iiz<iz+ioff+1; iiz++)
		    {
		      /* maintain periodic indices */
		      j = iiz;
		      if(iiz < 0)       j = NBINS+iiz;	       
		      if(iiz > NBINS-1) j = iiz-NBINS;
		      pixel_index = j + NBINS*iproc;
		      
		      rgrid =  (double)(j) * drbin;
		      
		      if (ilos[iproc] == 1)
			dr = fabs(rgrid-xx);
		      else if (ilos[iproc] == 2)
			dr = fabs(rgrid-yy);
		      else if (ilos[iproc] == 3)
			dr = fabs(rgrid-zz);
		      
		      if (dr > box2) dr = box100 - dr;
		      
		      dist2 = dr2 + dr*dr; /* dx^2+dy^2+dz^2 */		 
		      
		      if (dist2 <= h4)
			{
			  q = sqrt(dist2 * hinv2); /* q=r/h */

			  kernel = (q <= 1.0) 
			    ? 0.318309886-0.238732414*q*q*(2.0-q) : 0.079577471*(2.0-q)*(2.0-q)*(2.0-q);
			  
			  kernel *= hinv3;    /* kernel/h^3 */
			  
			  kernel  = dscale * kernel * P[i].Mass; /* g cm^-3 */
			  velker  = vscale * kernel * P[i].Vel[ilos[iproc]-1]; /* g cm^-3 * km s^-1 */
			  tempker = escale * kernel * P[i].U; /* g cm^-3 * K */
			  
			  rhoker_H[pixel_index]    += kernel * XH;		 
			  
			  rhoker_H1[pixel_index]   += kernel  * XH * P[i].NH0;
			  velker_H1[pixel_index]   += velker  * XH * P[i].NH0;
			  tempker_H1[pixel_index]  += tempker * XH * P[i].NH0;
			  
			} /* r < 2h */
		    } /* contributing vertices */
		} /* dx^2+dy^2 < (2h)^2 */
	    } /* dx < 2h */
	} /* particles in LOS */
      
      
      
      pcdone = 100.0*(float)iproc/(float)NUMLOS;
      if(pcdone >= pccount)
	{	  
	  printf("%3.2f%%\n",pcdone);
	  pccount += 10.0;
	}
      
    } /* LOS */
  
  
  return;
}

/********************************************************************************/

void get_los_positions()
{
  
  double Hz,drbin,vmax,dvbin,rscale;
  
  int i;
  
  Hz = 100.0*h100*sqrt(omega0/(atime*atime*atime)+omegaL); /* km s^-1 Mpc^-1 */
  
  drbin  = box100/(double)NBINS; /* comoving kpc/h */
  rscale = (KPC*atime)/h100;     /* convert length to cm */
  vmax   = box100*Hz*rscale/MPC; /* box size */
  dvbin  = vmax/(double)NBINS;   /* bin size (kms^-1) */

  for(i=0;i<NBINS-1;i++)
    {
      posaxis[i+1] = posaxis[i] + drbin; /* comoving kpc/h */
      velaxis[i+1] = velaxis[i] + dvbin; /* km s^-1 */
    }
  
   
#ifdef LOSTAB /* External ASCII file with line-of-sight positions */
  
  char fname[400];
  FILE *input;
  
  sprintf(fname,"%s",LOSFILE);
  if(!(input = fopen(fname, "r")))
    {
      printf("-DLOSTAB: loading external line-of-sight table\n");
      printf("Cannot read file %s\n\n",fname);
      exit(0);
    }
  printf("-DLOSTAB: loading external line-of-sight table\n");
  printf("Loaded %s\n\n",fname);

  for(i=0;i<NUMLOS;i++)                                            
    {
      fscanf(input, "%d %lf %lf %lf",&ilos[i],&xlos[i],&ylos[i],&zlos[i]);                                       
      xlos[i] *= 1.0e3;
      ylos[i] *= 1.0e3;
      zlos[i] *= 1.0e3;										       
    }
  
  fclose(input);

#else /* Random sight-lines or PRACE grid */
  
#ifdef PRACE_GRID
  
  /* Construct a regular grid for los coordinates.  This is the
     format used in the PRACE on-the-fly extraction */
  
  int ngrid1,ngrid2,ngrid3,ngrid_tot,iproc,row,col;
  double dgrid1,dgrid2,dgrid3;
  
  ngrid1 = 50;
  ngrid2 = 40;
  ngrid3 = 30;
  
  ngrid_tot = ngrid1*ngrid1 + ngrid2*ngrid2 + ngrid3*ngrid3;
  if(ngrid_tot != NUMLOS)
    {
      printf("\n\nError! When using PRACE_GRID, NUMLOS must be %d\n",ngrid_tot);
      printf("NUMLOS (see parameters.h) is currently %d\n\n",NUMLOS);
      exit(0);
    }
  
  dgrid1 = box100/(double)ngrid1;
  dgrid2 = box100/(double)ngrid2;
  dgrid3 = box100/(double)ngrid3;
  
  /* project along x-axis */
  for (row=0; row<ngrid1; row++)
    {
      for(col=0; col<ngrid1; col++)
	{
	  iproc = row*ngrid1 + col;
	  
	  ilos[iproc] = 1;
	  xlos[iproc] = 0.5*box100;
	  ylos[iproc] = (col + 0.131213)*dgrid1; /* Use random offset <0.5 to avoid ICs artifacts */
	  zlos[iproc] = (row + 0.131213)*dgrid1;
	}
    }
  /* project along y-axis */
  for (row=0; row<ngrid2; row++)
    {
      for(col=0; col<ngrid2; col++)
	{
          iproc = ngrid1*ngrid1 + row*ngrid2 + col;

          ilos[iproc] = 2;
          ylos[iproc] = 0.5*box100;
          xlos[iproc] = (col + 0.241008)*dgrid2;
          zlos[iproc] = (row + 0.241008)*dgrid2;
        }
    }
  /* project along z-axis */
  for (row=0; row<ngrid3; row++)
    {
      for(col=0; col<ngrid3; col++)
        {
          iproc = ngrid1*ngrid1 + ngrid2*ngrid2 + row*ngrid3 + col;

          ilos[iproc] = 3;
          zlos[iproc] = 0.5*box100;
          xlos[iproc] = (col + 0.170482)*dgrid3;
          ylos[iproc] = (row + 0.170482)*dgrid3;
        }
    }

#else
  
  /* Random lines-of-sight */
  srand48(241008); 
  
  for(i=0;i<NUMLOS;i++)  
    {
      do	
	ilos[i] = (int)(drand48()*4);
      while (ilos[i] == 0 || ilos[i] == 4); 
      
      xlos[i] = drand48()*box100;
      ylos[i] = drand48()*box100;
      zlos[i] = drand48()*box100;
    }
#endif
  
#endif
  
}


/********************************************************************************/

void normalise_fields()
{
  double H0,rhoc,critH;
  
  int i=0, iproc, pixel_index;
  
  H0    = 1.0e7/MPC; /* 100 km s^-1 Mpc^-1 in cgs */ 
  rhoc  = 3.0*(H0*h100)*(H0*h100)/(8.0*PI*GRAVITY); /* g cm^-3 */
  critH = rhoc*omegab*XH/(atime*atime*atime); /* g cm^-3*/

  /* Normalise quantities */
  for(iproc=0; iproc<NUMLOS; iproc++)
    {
      for(i = 0; i<NBINS; i++)
	{
	  pixel_index = i + NBINS*iproc;
	  
	  velker_H1[pixel_index]  /= rhoker_H1[pixel_index]; /* HI weighted km s^-1 */ 
	  tempker_H1[pixel_index] /= rhoker_H1[pixel_index]; /* HI weighted K */
	  rhoker_H1[pixel_index]  /= rhoker_H[pixel_index];  /* HI/H */
	  rhoker_H[pixel_index]   /= critH;
	}
    }
}

/********************************************************************************/


void compute_absorption()
{
  double T0,T1,T2;
  double rscale,escale,drbin,H0,rhoc,critH;
  double pcdone,pccount=0.0;
  double profile_H1,vdiff_H1;
  double sigma_Lya_H1,k1_H1,k2_H1,k3_H1;
  double u_H1[NBINS],binv_H1[NBINS],b2inv_H1[NBINS];
  double aa_H1[NBINS],tau_delta_H1[NBINS];
  
  int i,j,iproc,convol_index,pixel_index;
  
  rscale = (KPC*atime)/h100;     /* comoving kpc/h to cm */
  escale = 1.0e10;               /* (km s^-1)^2 to (cm s^-1)^2 */
  drbin  = box100/(double)NBINS; /* comoving kpc/h */
  H0     = 1.0e7/MPC; /* 100 km s^-1 Mpc^-1 in cgs */ 
  rhoc   = 3.0*(H0*h100)*(H0*h100)/(8.0*PI*GRAVITY); /* g cm^-3 */
  critH  = rhoc*omegab*XH/(atime*atime*atime); /* g cm^-3*/
  
  /* HI Lyman-alpha */
  sigma_Lya_H1 = sqrt(3.0*PI*SIGMA_T/8.0)*LAMBDA_LYA_H1 *FOSC_LYA; /* cm^2 */
  k1_H1        = 2.0*BOLTZMANN/(HMASS*AMU);
  k2_H1        = GAMMA_LYA_H1*LAMBDA_LYA_H1/(4.0*PI);
  k3_H1        = sigma_Lya_H1*C*rscale*drbin*critH/(sqrt(PI)*HMASS*AMU);

  for(iproc=0; iproc<NUMLOS; iproc++)
    {
      for(j=0; j<NBINS; j++)
  	{
  	  convol_index =  j + NBINS*iproc;
	  
  	  /* HI Lyman-alpha */
  	  u_H1[j]         = velaxis[j] + velker_H1[convol_index]; /* km s^-1 */
  	  binv_H1[j]      = 1.0/sqrt(k1_H1*tempker_H1[convol_index]); /* cm s^-1 */
  	  b2inv_H1[j]     = escale*binv_H1[j]*binv_H1[j];  /* (km s^-1)^-2 */
  	  aa_H1[j]        = k2_H1*binv_H1[j];
  	  tau_delta_H1[j] = k3_H1*rhoker_H[convol_index]*rhoker_H1[convol_index]*binv_H1[j];
  	}
      
      /* Voigt profile: Tepper-Garcia, 2006, MNRAS, 369, 2025 */
      for(i=0; i<NBINS; i++)
  	{
  	  pixel_index =  i + NBINS*iproc;
	  
  	  for(j=0; j<NBINS; j++)
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
      
       
      pcdone = 100.0*(float)iproc/(float)NUMLOS;
      if(pcdone >= pccount)
	{	  
	  printf("%3.2f%%\n",pcdone);
	  pccount += 10.0;
	}
           
    }
}
  

/********************************************************************************/

void write_spectra()
{
  double ztime[1],om_out[1],ol_out[1],ob_out[1],h_out[1],box_out[1],xh_out[1];
  int nbins_spec[1],nlos_spec[1];
  
  char fname[400];

  FILE *output;
  
  
  ztime[0]   = 1.0/atime - 1.0;
  om_out[0]  = omega0;
  ol_out[0]  = omegaL;
  ob_out[0]  = omegab;
  h_out[0]   = h100;
  xh_out[0]  = XH;
  box_out[0] = box100;
  
  nbins_spec[0] = (int)NBINS;
  nlos_spec[0]  = (int)NUMLOS;
  
  /* Saves the output in the format used for the on-the-fly data, so
     that the LOS and optical depth data is split over two files  */
#ifdef LOSEXTRACT
  
  sprintf(fname, "los%d_n%d_z%.3f.dat",NBINS,NUMLOS,ztime[0]);
  if(!(output = fopen(fname,"wb")))
    {
      printf("\nCannot write file:\n");
      printf("%s\n", fname);
      exit(0);
    }

  /* Header */
  fwrite(ztime,sizeof(double),1,output);
  fwrite(om_out,sizeof(double),1,output);
  fwrite(ol_out,sizeof(double),1,output);
  fwrite(ob_out,sizeof(double),1,output);
  fwrite(h_out,sizeof(double),1,output);
  fwrite(box_out,sizeof(double),1,output);
  fwrite(xh_out,sizeof(double),1,output);
  fwrite(nbins_spec,sizeof(int),1,output);
  fwrite(nlos_spec,sizeof(int),1,output);
  
  /* Positions */
  fwrite(ilos,sizeof(int),NUMLOS,output);
  fwrite(xlos,sizeof(double),NUMLOS,output);
  fwrite(ylos,sizeof(double),NUMLOS,output);
  fwrite(zlos,sizeof(double),NUMLOS,output);
  fwrite(posaxis,sizeof(double),NBINS,output);           /* pixel positions, comoving kpc/h */
  fwrite(velaxis,sizeof(double),NBINS,output);           /* pixel positions, km s^-1 */
  
  fwrite(rhoker_H,sizeof(double),NBINS*NUMLOS,output);   /* gas overdensity  */
  
  /* HI Lyman-alpha */
  fwrite(rhoker_H1,sizeof(double),NBINS*NUMLOS,output);  /* n_HI/n_H */
  fwrite(tempker_H1,sizeof(double),NBINS*NUMLOS,output); /* T [K], HI weighted */
  fwrite(velker_H1,sizeof(double),NBINS*NUMLOS,output);  /* v_pec [km s^-1], HI weighted */
  
  fclose(output);
  
  sprintf(fname, "tau%d_n%d_z%.3f.dat",NBINS,NUMLOS,ztime[0]);
  if(!(output = fopen(fname,"wb")))
    {
      printf("\nCannot write file:\n");
      printf("%s\n", fname);
      exit(0);
    }
  fwrite(tau_H1,sizeof(double),NBINS*NUMLOS,output);
  fclose(output);
  

#else


  sprintf(fname, "spec%d_n%d_z%.3f.dat",NBINS,NUMLOS,ztime[0]);
  if(!(output = fopen(fname,"wb")))
    {
      printf("\nCannot write file:\n");
      printf("%s\n", fname);
      exit(0);
    }

  /* Header */
  fwrite(ztime,sizeof(double),1,output);
  fwrite(om_out,sizeof(double),1,output);
  fwrite(ol_out,sizeof(double),1,output);
  fwrite(ob_out,sizeof(double),1,output);
  fwrite(h_out,sizeof(double),1,output);
  fwrite(box_out,sizeof(double),1,output);
  fwrite(xh_out,sizeof(double),1,output);
  fwrite(nbins_spec,sizeof(int),1,output);
  fwrite(nlos_spec,sizeof(int),1,output);
  
  /* Positions */
  fwrite(ilos,sizeof(int),NUMLOS,output);
  fwrite(xlos,sizeof(double),NUMLOS,output);
  fwrite(ylos,sizeof(double),NUMLOS,output);
  fwrite(zlos,sizeof(double),NUMLOS,output);
  fwrite(posaxis,sizeof(double),NBINS,output);           /* pixel positions, comoving kpc/h */
  fwrite(velaxis,sizeof(double),NBINS,output);           /* pixel positions, km s^-1 */
  
  fwrite(rhoker_H,sizeof(double),NBINS*NUMLOS,output);   /* gas overdensity  */
  
  /* HI Lyman-alpha */
  fwrite(rhoker_H1,sizeof(double),NBINS*NUMLOS,output);  /* n_HI/n_H */
  fwrite(tempker_H1,sizeof(double),NBINS*NUMLOS,output); /* T [K], HI weighted */
  fwrite(velker_H1,sizeof(double),NBINS*NUMLOS,output);  /* v_pec [km s^-1], HI weighted */
  fwrite(tau_H1,sizeof(double),NBINS*NUMLOS,output);     /* HI optical depth */
  
  fclose(output);

#endif
}
  

/********************************************************************************/

void InitLOSMemory()
{  
  
  ilos      = (int *)calloc(NUMLOS, sizeof(int));
  if(NULL==ilos){free(ilos);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}

  xlos      = (double *)calloc(NUMLOS, sizeof(double));
  if(NULL==xlos){free(xlos);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}

  ylos      = (double *)calloc(NUMLOS, sizeof(double));
  if(NULL==ylos){free(ylos);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}

  zlos      = (double *)calloc(NUMLOS, sizeof(double));
  if(NULL==zlos){free(zlos);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}

  posaxis      = (double *)calloc(NBINS, sizeof(double));
  if(NULL==posaxis){free(posaxis);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}

  velaxis      = (double *)calloc(NBINS, sizeof(double));
  if(NULL==posaxis){free(posaxis);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}
  
  rhoker_H     = (double *)calloc(NUMLOS*NBINS, sizeof(double));
  if(NULL==rhoker_H){free(rhoker_H);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}
  
  rhoker_H1    = (double *)calloc(NUMLOS*NBINS, sizeof(double));
  if(NULL==rhoker_H1){free(rhoker_H1);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}

  velker_H1    = (double *)calloc(NUMLOS*NBINS, sizeof(double));
  if(NULL==velker_H1){free(velker_H1);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}

  tempker_H1   = (double *)calloc(NUMLOS*NBINS, sizeof(double));
  if(NULL==tempker_H1){free(tempker_H1);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}

  tau_H1       = (double *)calloc(NUMLOS*NBINS, sizeof(double));
  if(NULL==tau_H1){free(tau_H1);printf("Memory allocation failed in extract_spectra.c\n");exit(0);}

}


/********************************************************************************/
