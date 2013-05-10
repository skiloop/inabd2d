# compiler
CC=cc#icc#gcc

# linker
LD=ld

# source path
SRC=./src

VPATH= $(SRC)
#
# origin CFLAGS
CFLAGS=-O4 -Wall#-Wunused-but-set-variable# run mode
#CFLAGS=-g -Wall -DDEBUG#-Wunused-but-set-variable #debug mode

# origin link options
LIB=-lm
LIB+=-O4

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

# add -MMD
CFLAGS+=-MMD

OBJS=abd2d.o  \
commonData.o \
connectingInterface.o \
dataSaving.o \
dataType.o \
density.o \
fdtd.o \
initial.o \
breakdownFormula.o \
matlabSimulation.o \
pml.o 
.PHONY: all clean
all:abd2d
abd2d:$(OBJS)
	$(CC) -o abd2d $^ $(LIB)
$(OBJS):%.o:%.c common.h
	$(CC) $(CFLAGS) -c $<
clean:
	-rm -f *.o abd2d *.d


