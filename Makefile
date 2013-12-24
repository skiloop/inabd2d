
include makefile.in

EXCUTABLE:=abd2d
SRC_DIR:=./src
TEST_SRC_DIR:=./test
VPATH = $(SRC_DIR):$(TEST_SRC_DIR)
SOURCES=$(shell find $(SRC_DIR) -name "*.c")
OBJS:=$(patsubst $(SRC_DIR)/%.c,%.o,$(SOURCES))
DEPS:=$(patsubst $(SRC_DIR)/%.c,%.d,$(SOURCES))
INCLUDE+=-I$(SRC_DIR)
#We don't need to clean up when we're making these targets
NODEPS:=clean tags svn
TEST:=sourceTest
TARGET_OBJS:=abd2d.o sourceTest.o
NORMAL_OBJS:=$(filter-out $(TARGET_OBJS),$(OBJS))
PROJECTS=$(TEST) $(EXCUTABLE) #3DFormulaTransforming.pdf
.PHONY:all clean test objs veryclean rebuild deps
	
all:$(PROJECTS) 

deps:$(DEPS)

objs:$(OBJS)

test:$(TEST)

sourceTest:sourceTest.o $(NORMAL_OBJS)
	echo "make " $@
	$(CC) -o $@ $^ $(LIB)

#This is the rule for creating the dependency files
%.d:$(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -MM -MT $*.o $< -MF $@

%.o:%.c %.d
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $< 

#Don't create dependencies when we're cleaning, for instance
ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
#Chances are, these files don't exist.  GMake will create them and
#clean up automatically afterwards
-include $(DEPS)
endif

$(EXCUTABLE):$(OBJS)
	$(CC) -o $@ $(OBJS) $(LIB)
clean:
	-rm -f $(DEPS) $(OBJS) $(PROJECTS) *.aux *.log
veryclean:clean

rebuild:veryclean all

