# LIBRARY SETTINGS - SET AS NECESSARY

INC_BOOST = /home/dilthey/PnP/libs/boost_1_52_0/
LIB_BOOST = /home/dilthey/PnP/libs/boost_1_52_0/stage/lib/
INCS = -I$(INC_BOOST) 
LIBS = -L$(LIB_BOOST) -lboost_system -lboost_filesystem


MKDIR_P = mkdir -p

.PHONY: directories
	
# END LIBRARY SETTINGS

#
# object and binary dirs
#

DIR_OBJ = ../obj
DIR_BIN = ../bin

CXX    = g++
COPTS  = -ggdb -O2 -fopenmp -std=gnu++0x -fstack-protector-all
CFLAGS = 
COMPILE = $(CXX) $(INCS) $(CFLAGS) $(COPTS)

VPATH = Data:Graph:NextGen:hash:hash/deBruijn:hash/sequence:GraphAligner:GraphAlignerUnique:readFilter
        
OBJS = \
        $(DIR_OBJ)/Utilities.o \
        $(DIR_OBJ)/DeBruijnGraph.o \
        $(DIR_OBJ)/DeBruijnElement.o \
        $(DIR_OBJ)/basic.o \
        $(DIR_OBJ)/binarykMer.o \
        $(DIR_OBJ)/Hsh.o \
        $(DIR_OBJ)/Validator.o \
#
# list executable file names
#
EXECS = kMerificator

OUT_DIR = ../obj ../bin

directories: ${OUT_DIR}


#
# compile and link
#
default:
	@echo
	@echo " to build:"
	@echo "    make all"
	@echo
	@echo " to clean:"
	@echo "    make clean"
	@echo "    make realclean"
	@echo

all: directories $(EXECS)

$(EXECS): $(OBJS)
	$(foreach EX, $(EXECS), $(COMPILE) $(EX).cpp -c -o $(DIR_OBJ)/$(EX).o;)
	$(foreach EX, $(EXECS), $(COMPILE) $(OBJS) $(DIR_OBJ)/$(EX).o -o $(DIR_BIN)/$(EX) $(LIBS);)

$(DIR_OBJ)/%.o: %.cpp %.h
	$(COMPILE) $< -c -o $@


#
# odds and ends
#
clean:
	/bin/rm $(DIR_OBJ)/*

realclean: clean
	/bin/rm $(DIR_BIN)/*

${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}

