
#######################################################################
#  Look at end of file for a brief guide to any compile-time options. #
#######################################################################

#--------------------------------------- Spectra extraction options

OPTS +=-DLOSEXTRACT	
OPTS +=-DPRACE_GRID	
#OPTS +=-DLOSTAB     

#--------------------------------------- Select system
#SYSTYPE="Darwin"
#SYSTYPE="pppzjsb"
#SYSTYPE="dator"
#SYSTYPE="brahan"
SYSTYPE="ulysses"

ifeq ($(SYSTYPE),"ulysses")
CC = gcc
OPTIMIZE = -O3
CFLAGS =  $(OPTIMIZE) -Wall 
HDF5INCL = -I/scratch/rmurgia/miniconda2/include/
HDF5LIB  = -L/scratch/rmurgia/miniconda2/lib/ -lhdf5 -lz
endif

ifeq ($(SYSTYPE),"Darwin")
CC = icc
OPTIMIZE = -O3
CFLAGS =  $(OPTIMIZE) -Wall -xHOST -ip -ipo
HDF5INCL = -I/usr/local/Cluster-Apps/hdf5/1.8.0/include/
HDF5LIB  = -L/usr/local/Cluster-Apps/hdf5/1.8.0/lib/ -lhdf5 -lz
endif

ifeq ($(SYSTYPE),"pppzjsb")
CC = gcc		
OPTIMIZE = -O3
CFLAGS = $(OPTIMIZE) -Wall
HDF5INCL = -I/home/ppzjsb/FILES/Libraries/hdf5-1.8.11-linux-x86_64-shared/include
HDF5LIB  = -L/home/ppzjsb/FILES/Libraries/hdf5-1.8.11-linux-x86_64-shared/lib/ -lhdf5 -lz
endif

ifeq ($(SYSTYPE),"dator")
CC = gcc		
OPTIMIZE = -O3
CFLAGS = $(OPTIMIZE) -Wall
HDF5INCL = -I/home/ppzjsb/Libraries/hdf5-1.8.11-linux-x86_64-shared/include
HDF5LIB  = -L/home/ppzjsb/Libraries/hdf5-1.8.11-linux-x86_64-shared/lib/ -lhdf5 -lsz 
endif

ifeq ($(SYSTYPE),"brahan")
CC = gcc		
OPTIMIZE = -O3
CFLAGS = $(OPTIMIZE) -Wall
HDF5INCL = -I/home/ppzjsb/Libraries/hdf5-1.8.11-linux-x86_64-shared/include
HDF5LIB  = -L/home/ppzjsb/Libraries/hdf5-1.8.11-linux-x86_64-shared/lib/ -lhdf5 -lsz 
endif


CFLAGS += $(OPTS) $(HDF5INCL)

EXEC = SpecExtract
OBJS = read_snapshot.o extract_spectra.o
INCL = parameters.h global_vars.h proto.h 

LIBS = -lm $(HDF5LIB)


$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)
$(OBJS): $(INCL)

.PHONY: clean
clean:
	rm -f $(OBJS) $(EXEC)

##############################################################################
#
#
# Spectra extraction options:
#
#	- LOSEXTRACT 	Saves the output in the format used for the on-the-fly
#			data, so that the LOS and optical depth data is split 
#			over two files 
#
#	- PRACE_GRID	Extract 5000 sight-lines on a regular grid, 50^2, 40^2,
#			30^2 along each axis. These are the positions used in 
#			the PRACE on-the-fly extraction 
#
#	- LOSTAB	pre-define line-of-sight positions in an external 
# 			table (see parameters.h)	
#
#
##############################################################################
