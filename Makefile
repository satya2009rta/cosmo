####################
# compilers
####################
CC        = g++
CXXFLAGS 	= -Wall -Wextra -std=c++2a -O3 -DNDEBUG -Wno-unknown-pragmas 
CFLAGS 	= -O3 -DNDEBUG

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

.PHONY: folder clean cosmo genMaze

TARGET = folder cosmo genMaze

build: $(TARGET)

folder:
	[ -d $(BUILD) ] || mkdir -p $(BUILD)

clean:
	rm -r -f  $(BUILD)/*

cosmo:
	$(CC) $(CXXFLAGS) $(LIBINC) $(SRC)/cosmo.cpp -o $(BUILD)/cosmo

genMaze:
	$(CC) $(CXXFLAGS) $(LIBINC) $(BSRC)/genMaze.cpp -o $(BUILD)/genMaze