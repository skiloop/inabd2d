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
breakdownParamters.o \
commonData.o \
connectingInterface.o \
dataSaving.o \
fdtd.o \
fieldCreation.o \
initCommonData.o \
init.o \
InonizationFormula.o \
matlabSimulation.o \
memoryRelease.o \
myStruct.o \
pml.o \
updateDensity.o
.PHONY: all clean
all:abd2d
abd2d:$(OBJS)
	$(CC) -o abd2d $^ $(LIB)
$(OBJS):%.o:%.c common.h
	$(CC) $(CFLAGS) -c $<
clean:
	-rm -f *.o abd2d *.d


