
/* Function prototypes */

/* read_snapshot.c */
int allocate_memory();
int load_snapshot_hdf5(int i);


/* extract_spectra.c */
void SPH_interpolation();
void InitLOSMemory();
void compute_absorption();
void normalise_fields();
void write_spectra();
void get_los_positions();
