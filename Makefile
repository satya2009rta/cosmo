####################
# compilers
####################
CC        = g++
C	  		= gcc
CXXFLAGS 	= -Wall -Wextra -std=c++2a -O3 -DNDEBUG -Wno-unknown-pragmas 
CFLAGS 	= -O3 -DNDEBUG


## OpenMP compiler and flags
# MPCC 	= g++
# MPFLAGS	= -fopenmp



# No OpenMP
MPCC 		= g++
MPFLAGS		=

####################
# project root
####################
PRJROOT		= .
LIBINC		= -I$(PRJROOT)/lib -isystem$(PRJROOT)
SRC			= $(PRJROOT)/src
BSRC		= $(PRJROOT)/benchmarking/src
BUILD		= $(PRJROOT)/build


#
# main executable files
#

.PHONY: folder clean compute negotiate pg2dpg genMaze

TARGET = folder compute negotiate pg2dpg genMaze

build: $(TARGET)

folder:
	[ -d $(BUILD) ] || mkdir -p $(BUILD)

clean:
	rm -r -f  $(BUILD)/*

compute:
	$(CC) $(CXXFLAGS) $(LIBINC) $(SRC)/compute.cpp -o $(BUILD)/compute

negotiate:
	$(MPCC) $(CXXFLAGS) $(LIBINC) $(MPFLAGS) $(SRC)/negotiate.cpp -o $(BUILD)/negotiate

pg2dpg:
	$(CC) $(CXXFLAGS) $(LIBINC) $(SRC)/pg2dpg.cpp -o $(BUILD)/pg2dpg

genMaze:
	$(CC) $(CXXFLAGS) $(LIBINC) $(BSRC)/genMaze.cpp -o $(BUILD)/genMaze