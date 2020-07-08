#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>


#include "global_vars.h"
#include "parameters.h"
#include "proto.h"



char *path,*snapbase;
int  snapshot;



/********************************************************************************/

int main(int argc, char **argv)
{
  int i;

  printf("\n");
  
  if(argc != 4)
    {
      fprintf(stderr, "\n\nIncorrect argument(s).  Specify:\n\n");
      fprintf(stderr, "<path>     (path to output files)\n");
      fprintf(stderr, "<snapbase> (basename for snapshots)\n");
      fprintf(stderr, "<num>      (number of snapshot)\n\n");
      fprintf(stderr, "Usage: ./SpecExtract <path> <snapbase> <num>\n");
      fprintf(stderr, "\n\n");
      exit(0);
    }
  
  path     = argv[1];
  snapbase = argv[2];
  snapshot = atoi(argv[3]);

#if defined(PRACE_GRID) && defined(LOSTAB)
  printf("The flags PRACE_GRID and LOSTAB cannot be raised simultaenously\n");
  printf("Adjust the parameter settings in Makefile\n\n");
  exit(0);
#endif
  
  /* Load multiple files and project */
  for(i=0; i<FILENUMBER; i++) 
    { 
      load_snapshot_hdf5(i);
      
      if(i == 0)
	{
	  InitLOSMemory();
	  get_los_positions();
	  printf("Assuming hydrogen fraction by mass of %f\n\n",XH);
	}
      
      printf("Extracting LOS from sub-file %d...\n",i);
      SPH_interpolation();
      printf("Done.\n\n");
      
      free(P);
    }
  
  printf("Omegab0 = %f\n\n",omegab);
  printf("Normalising...\n");
  normalise_fields();
  printf("Done.\n\n");
  
  printf("Computing optical depths...\n");
  compute_absorption();
  printf("Done.\n\n");
  
  printf("Saving spectra...\n");
  write_spectra();
  printf("Done.\n\n");
  
  return(0);
}


/********************************************************************************/

/* loads the particle data from HDF5 snapshot files */
int load_snapshot_hdf5(int i)
{
  FILE *fd;
  char   buf[400];
  int    k,n,flag_cooling;
  hid_t hdf5_file, hdf5_grp, hdf5_attribute, hdf5_dataset, dapl_id=0;
  float *dummy3_float, *dummy_float;

  static double fgas=0.0, fdark=0.0, fstars=0.0, ftot=0.0;
  
  
  if(FILENUMBER>1)
    sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d.hdf5", path,snapshot,snapbase,snapshot,i);
  else
    sprintf(buf, "%s/%s_%03d.hdf5", path, snapbase, snapshot);
  
  if(!(fd=fopen(buf,"r")))
    {
      printf("can't open file `%s`\n\n",buf);
      exit(0);
    }
  
  /* Open HDF5 file */
  hdf5_file = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
  
  /* Read the header data*/
  hdf5_grp = H5Gopen(hdf5_file, "/Header", H5P_DEFAULT);
  
  hdf5_attribute = H5Aopen_name(hdf5_grp, "NumPart_ThisFile");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, NumPart);
  H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_grp, "NumPart_Total");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, NumPartTot);
  H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_grp, "NumPart_Total_HighWord");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, NumPartTotHigh);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_grp, "MassTable");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, mass);
  H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_grp, "Time");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &atime);
  H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_grp, "Redshift");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &redshift);
  H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_grp, "BoxSize");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &box100);
  H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_grp, "Omega0");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &omega0);
  H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_grp, "OmegaLambda");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &omegaL);
  H5Aclose(hdf5_attribute);
      
  hdf5_attribute = H5Aopen_name(hdf5_grp, "HubbleParam");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &h100);
  H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_grp, "Flag_Cooling");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, &flag_cooling);
  H5Aclose(hdf5_attribute);
  
  H5Gclose(hdf5_grp);
  
  
  if(i == 0)
    {
      for(k=0; k<6; k++)	
	{
	  printf("Npart[%d]        = %d\n",k,NumPart[k]);
	  printf("NpartTot[%d]     = %u\n",k,NumPartTot[k]);
	  printf("NpartTotHigh[%d] = %u\n",k,NumPartTotHigh[k]);
	  printf("Massarr[%d]      = %e\n\n",k,mass[k]);
	}
      
      printf("Redshift = %f\n",redshift);
      printf("a(t)     = %f\n",atime);
      printf("Lambda   = %f\n",omegaL);
      printf("Matter   = %f\n",omega0);
      printf("Hubble   = %f\n",h100);
      printf("Box      = %f\n\n",box100);
    }


    
  /* Allocate memory for particles in this sub-file */
  Ntype    = NumPart[0];
  NtypeTot = (long long)NumPartTot[0] + ((long long)NumPartTotHigh[0] << 32);
  allocate_memory();
  
  /* Gas particles */
  dummy3_float = (float *)calloc(Ntype*3, sizeof(float));
  if(NULL==dummy3_float){free(dummy3_float);printf("Memory allocation failed in read_snapshot.c\n");exit(0);}
  
  hdf5_grp = H5Gopen(hdf5_file, "/PartType0", H5P_DEFAULT);
      
  hdf5_dataset = H5Dopen(hdf5_grp, "Coordinates", dapl_id);
  H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dummy3_float);
  H5Dclose(hdf5_dataset);
  for(n=0; n<Ntype; n++)
    {  
      P[n].Pos[0] = dummy3_float[3*n];
      P[n].Pos[1] = dummy3_float[3*n+1];
      P[n].Pos[2] = dummy3_float[3*n+2];
    }
  if(i==0) 
    { 
      printf("P[0].Pos[0] = %f\n",  P[0].Pos[0]);
      printf("P[0].Pos[1] = %f\n",  P[0].Pos[1]);
      printf("P[0].Pos[2] = %f\n",  P[0].Pos[2]);
    }
  

  hdf5_dataset = H5Dopen(hdf5_grp, "Velocities", dapl_id);
  H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dummy3_float);
  H5Dclose(hdf5_dataset);
  for(n=0; n<Ntype; n++)
    {  
      P[n].Vel[0] = dummy3_float[3*n];
      P[n].Vel[1] = dummy3_float[3*n+1];
      P[n].Vel[2] = dummy3_float[3*n+2];
    }
  if(i==0) 
    { 
      printf("P[0].Vel[0] = %f\n", P[0].Vel[0]);
      printf("P[0].Vel[1] = %f\n", P[0].Vel[1]);
      printf("P[0].Vel[2] = %f\n", P[0].Vel[2]); 
    }
  
  free(dummy3_float);
  
  
  dummy_float = (float *)calloc(Ntype, sizeof(float));
  if(NULL==dummy_float){free(dummy_float);printf("Memory allocation failed in read_snapshot.c\n");exit(0);}
  
  hdf5_dataset = H5Dopen(hdf5_grp, "Masses", dapl_id);
  H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dummy_float);
  H5Dclose(hdf5_dataset);
  for(n=0; n<Ntype; n++)
    { 
      P[n].Mass = dummy_float[n];
      fgas += dummy_float[n];
    }
  if(i==0) printf("P[0].Mass   = %e\n", P[0].Mass);
  
  
  
  hdf5_dataset = H5Dopen(hdf5_grp, "InternalEnergy", dapl_id);
  H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dummy_float);
  H5Dclose(hdf5_dataset);
  for(n=0; n<Ntype; n++)
    P[n].U = dummy_float[n];
  if(i==0) printf("P[0].U      = %f\n", P[0].U);
      
      
  
  if(flag_cooling)
    {
      hdf5_dataset = H5Dopen(hdf5_grp, "ElectronAbundance", dapl_id);
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dummy_float);
      H5Dclose(hdf5_dataset);
      for(n=0; n<Ntype; n++)
	P[n].Ne = dummy_float[n];
      if(i==0) printf("P[0].Ne     = %e\n", P[0].Ne);
	  
      hdf5_dataset = H5Dopen(hdf5_grp, "NeutralHydrogenAbundance", dapl_id);
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dummy_float);
      H5Dclose(hdf5_dataset);
      for(n=0; n<Ntype; n++)
	P[n].NH0 = dummy_float[n];
      if(i==0) printf("P[0].NH0    = %e\n", P[0].NH0);
    }
  
  hdf5_dataset = H5Dopen(hdf5_grp, "SmoothingLength", dapl_id);
  H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dummy_float);
  H5Dclose(hdf5_dataset);
  for(n=0; n<Ntype; n++)
    P[n].h = dummy_float[n];
  if(i==0) printf("P[0].h      = %f\n\n", P[0].h);
  H5Gclose(hdf5_grp);

      
  /* Star particles */
  if (NumPart[4] > 0)
    {
      hdf5_grp = H5Gopen(hdf5_file, "/PartType4", H5P_DEFAULT);
      hdf5_dataset = H5Dopen(hdf5_grp, "Masses", dapl_id);
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dummy_float);
      H5Dclose(hdf5_dataset);
      H5Gclose(hdf5_grp);
      
      for(n=0; n<NumPart[4];n++)
	fstars += dummy_float[n];
    }
  
  H5Fclose(hdf5_file);
  free(dummy_float);
  
  fdark  += NumPart[1]*mass[1];
  
  ftot = fdark+fgas+fstars;
  
  if(i == FILENUMBER-1)
    omegab = (fgas+fstars)/ftot * omega0;
    
  return(0);
}


/********************************************************************************/

/* this routine allocates the memory for the particle data */
int allocate_memory()
{ 
  if(!(P=calloc(Ntype,sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n\n");
      exit(0);
    }
  
  return(0);
}

/********************************************************************************/
