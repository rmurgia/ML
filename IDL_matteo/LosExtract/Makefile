
#######################################################################
#  Look at end of file for a brief guide to any compile-time options. #
#######################################################################

#--------------------------------------- Select system
#SYSTYPE="Darwin"
#SYSTYPE="pppzjsb"
SYSTYPE="brahan"


ifeq ($(SYSTYPE),"Darwin")
CC = icc
OPTIMIZE = -O3
CFLAGS =  $(OPTIMIZE) -Wall -xHOST -ip -ipo
CFLAGS += $(OPTS)	
endif

ifeq ($(SYSTYPE),"pppzjsb")
CC = gcc	
OPTIMIZE = -O3
CFLAGS =  $(OPTIMIZE) -Wall
CFLAGS += $(OPTS)	
endif

ifeq ($(SYSTYPE),"brahan")
CC = gcc	
OPTIMIZE = -O3
CFLAGS =  $(OPTIMIZE) -Wall
CFLAGS += $(OPTS)	
endif

EXEC = LosExtract
OBJS = extract_los.o
INCL = proto.h global_vars.h
LIBS = -lm


$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)
$(OBJS): $(INCL)

.PHONY: clean
clean:
	rm -f $(OBJS) $(EXEC)
