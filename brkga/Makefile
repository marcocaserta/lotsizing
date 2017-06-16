#####################################################
# Author : Marco Caserta
# E-mail : marco dot caserta at ie dot edu
# Date   : 30.01.15
#
# This makefile is used to compile the 
# sources (c++) that implements the
# algorithm for the MMKP
#
#####################################################
SRCDIR        = src
BINDIR        = bin
OBJDIR        = obj

# ---------------------------------------------------------------------
# ILOG CPLEX stuff
#SYSTEM     = x86-64_debian4.0_4.1
SYSTEM     = x86-64_linux


LIBFORMAT  = static_pic
# ---------------------------------------------------------------------
# CPLEXDIR      = /home/ubuntu/opt/ibm/ILOG/cplex1261/cplex
# CONCERTDIR    = /home/ubuntu/opt/ibm/ILOG/cplex1261/concert
# CONCERTDIR    = /home/marco/opt/ibm/ILOG/cplex127/concert
# CPLEXDIR      = /home/marco/opt/ibm/ILOG/cplex127/cplex
CPLEXDIR      = /home/marco/opt/ibm/cplex1261/cplex
CONCERTDIR    = /home/marco/opt/ibm/cplex1261/concert
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
# ---------------------------------------------------------------------
# END ILOG CPLEX stuff

CCOPT     = -m64 -O -fPIC -fexceptions  -DIL_STD  -DOPTIMAL 
#-DNDEBUG
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread
CCFLAGS   = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 
#CTIME     = -pg -fprofile-arcs

# set default name for the executable
EXEC ?= mmkp

# set compiler
CC = g++

# set debug options
DEBUG = -ggdb

#set optimization level
OPTLEVEL = -O -DEBIAN_BUILDARCH=pentium

#set flags and warnings
FLAGS = $(DEBUG) $(OPTLEVEL) -DIL_STD -fexceptions -fPIC -pipe -Wparentheses -Wreturn-type  -Wpointer-arith -Wwrite-strings -Wconversion -Wno-vla 

#-Wcast-qual
#-fexceptions
#-fPIC
#-DNDEBUG     # <<<---- Activate it to deactivate all the assert()
#-Wall
#-fomit-frame-pointer
#-Wimplicit

#default: $(OBJDIR)/SampleDecoder.o $(OBJDIR)/lagrange.o $(OBJDIR)/cplex.o $(OBJDIR)/mmkp.o $(OBJDIR)/options.o $(OBJDIR)/timer.o 
#	$(CC) -g $(CCFLAGS) $(RINCLUDE) $(CTIME) $(OBJDIR)/SampleDecoder.o $(OBJDIR)/options.o $(OBJDIR)/lagrange.o $(OBJDIR)/cplex.o $(OBJDIR)/timer.o $(OBJDIR)/mmkp.o -o $(BINDIR)/$(EXEC)# $(CCLNFLAGS)
default: $(OBJDIR)/lagrange.o $(OBJDIR)/cplex.o $(OBJDIR)/mmkp.o $(OBJDIR)/options.o $(OBJDIR)/timer.o $(OBJDIR)/robust.o 
	$(CC) -g $(CCFLAGS) $(RINCLUDE) $(CTIME) $(OBJDIR)/options.o $(OBJDIR)/lagrange.o $(OBJDIR)/cplex.o $(OBJDIR)/timer.o $(OBJDIR)/mmkp.o $(OBJDIR)/robust.o -o $(BINDIR)/$(EXEC) $(CCLNFLAGS)

$(OBJDIR)/mmkp.o: $(SRCDIR)/mmkp.cpp 
	$(CC) -c $(CCFLAGS) $(CTIME) $(SRCDIR)/mmkp.cpp -o $(OBJDIR)/mmkp.o
$(OBJDIR)/lagrange.o: $(SRCDIR)/lagrange.cpp
	$(CC) -c $(CCFLAGS) $(CTIME) $(SRCDIR)/lagrange.cpp -o $(OBJDIR)/lagrange.o 
$(OBJDIR)/cplex.o: $(SRCDIR)/cplex.cpp
	$(CC) -c $(CCFLAGS) $(CTIME) $(SRCDIR)/cplex.cpp -o $(OBJDIR)/cplex.o 
$(OBJDIR)/robust.o: $(SRCDIR)/robust.cpp
	$(CC) -c $(CCFLAGS) $(CTIME) $(SRCDIR)/robust.cpp -o $(OBJDIR)/robust.o 
$(OBJDIR)/options.o: $(SRCDIR)/options.cpp
	$(CC) -c -g $(FLAGS) $(CTIME) $(SRCDIR)/options.cpp -o $(OBJDIR)/options.o
$(OBJDIR)/timer.o: $(SRCDIR)/timer.cpp
	$(CC) -c -g $(FLAGS) $(CTIME) $(SRCDIR)/timer.cpp -o $(OBJDIR)/timer.o

brkga: $(OBJDIR)/SampleDecoder.o $(OBJDIR)/brkgaMain.o
	$(CC) -g $(OBJDIR)/brkgaMain.o $(OBJDIR)/SampleDecoder.o -o $(BINDIR)/brkga
$(OBJDIR)/brkgaMain.o: $(SRCDIR)/brkgaMain.cpp
	$(CC) -c $(SRCDIR)/brkgaMain.cpp -o $(OBJDIR)/brkgaMain.o 
$(OBJDIR)/SampleDecoder.o: $(SRCDIR)/SampleDecoder.cpp
	$(CC) -c $(SRCDIR)/SampleDecoder.cpp -o $(OBJDIR)/SampleDecoder.o 

# clean backup files
clean:
	rm *.*~ 
	rm *.o

# create doxygen documentation using "doxygen.conf" file
# the documentation is put into the directory Doc	
doc: $(SRCDIR)/lad.cpp doxygen.conf
	doxygen doxygen.conf


#	
