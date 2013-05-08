# compiler
CC= gcc #icc#gcc

# linker
LD=ld

# current path
CURDIR=.

#
# origin CFLAGS
CFLAGS=-O3 -Wall #-Wunused-but-set-variable# run mode
#CFLAGS=-g -Wall -DDEBUG #-Wunused-but-set-variable #debug mode

# origin link options
LIB=-lm

# source path
SRC=$(CURDIR)/src

# matlab path 
MATPATH=/opt/Matlab/R2011a# for local server
#MATPATH=/opt2/Matlab/R2011a# for servers of college

# Matlab link path
MATLINK=-Wl,-rpath=$(MATPATH)/bin/glnxa64

# Matlab link option
MATLIB=$(MATLINK) -L$(MATPATH)/bin/glnxa64 -lmx -leng

# Matlab link path
MATINC=-I$(MATPATH)/extern/include

# link option
#LIB+=$(MATLIB)

# add Matlab simulation
#CFLAGS+=-DMATLAB_SIMULATION $(MATINC)# -g

all:abd2d
abd2d:abd2d-niu.o InonizationFormula.o
	$(CC) -o abd2d abd2d-niu.o InonizationFormula.o $(LIB)
abd2d-niu.o:$(SRC)/abd2d-niu.c $(SRC)/*.h
	$(CC) $(CFLAGS) -c $(SRC)/abd2d-niu.c
InonizationFormula.o:$(SRC)/InonizationFormula.c $(SRC)/*.h
	$(CC) $(CFLAGS) -c $(SRC)/InonizationFormula.c
clean:
	rm -f *.o abd2d


