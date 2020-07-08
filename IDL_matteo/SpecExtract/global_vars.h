

struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  
  float  Rho, U, h;
  float NH0, NHep, Ne;
} *P;


int  NumPart[6], Ntype;
unsigned int NumPartTot[6], NumPartTotHigh[6];
long long NtypeTot;
double  atime, redshift, omega0, omegaL, omegab, box100, h100, mass[6];



/* Numbers */
#define  PI    3.14159265358979323846
#define  GAMMA (5.0/3.0)

/* Physical constants (cgs units) */
/* See http://physics.nist.gov/cuu/Constants/index.html */ 
#define  GRAVITY     6.67384e-8
#define  BOLTZMANN   1.3806488e-16
#define  C           2.99792458e10
#define  AMU         1.66053886e-24 /* 1 a.m.u */
#define  MPC         3.08568025e24
#define  KPC         3.08568025e21
#define  SIGMA_T     6.652458734e-25 
#define  SOLAR_MASS  1.989e33

/* Atomic data (from VPFIT) */
#define  LAMBDA_LYA_H1  1215.6701e-8 /* cm */
#define  FOSC_LYA       0.416400
#define  HMASS          1.00794  /* Hydrogen mass in a.m.u. */
#define  GAMMA_LYA_H1   6.265e8  /* s^-1 */


#define XH 0.76  /* hydrogen fraction by mass */
