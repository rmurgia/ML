#include <stdio.h>

/*---------------------- Global parameters ---------------------------*/

/* Snapshot information */
#define FILENUMBER 256 /* Number of snapshot sub-files */

/* Spectra information */
#define NBINS  2048 /* number of pixels in each line-of-sight*/
#define NUMLOS 5000   /* number of lines-of-sight */


/* --------------------- OPTS +=-DLOSTAB -----------------------------*/

/* Path and name for external table with LOS co-ordinates */
/* The table format is: axis(=1,2 or 3), x, y, z [cMpc/h] */

#define LOSFILE  "./table_los_grid.txt"

/*--------------------------------------------------------------------*/
