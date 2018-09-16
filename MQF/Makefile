TARGETS=main
TESTFILES = tests/CountingTests.o tests/HighLevelFunctionsTests.o tests/IOTests.o tests/tagTests.o

ifdef D
	DEBUG=-g
	OPT=
else
	DEBUG=
	OPT=-O3
	#OPT= -O3 -no-prec-div -ansi-alias
endif

ifdef NH
	ARCH=
else
	ARCH=
endif

ifdef P
	PROFILE=-pg -no-pie # for bug in gprof.
endif

CXX = g++
CC = g++
LD= g++

CXXFLAGS =   -std=c++11 -fPIC -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I. -Wno-unused-result -Wno-strict-aliasing -Wno-unused-function -fopenmp

LDFLAGS = $(DEBUG) $(PROFILE) $(OPT)

#
# declaration of dependencies
#

all: $(TARGETS)



# dependencies between programs and .o files

main:	main.o	gqf.o	utils.o
	$(LD) $^ $(LDFLAGS) -o $@
# dependencies between .o files and .h files

libgqf: gqf.o utils.o
	$(LD) $^ $(LDFLAGS) --shared -o $@.so

test:  $(TESTFILES) gqf.c test.o utils.o
	$(LD) $(LDFLAGS) -DTEST -o mqf_test test.o utils.o $(TESTFILES) gqf.c

main.o: 								 									gqf.h

# dependencies between .o files and .cc (or .c) files


gqf.o: gqf.c gqf.h

#
# generic build rules
#







%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -c -o $@

%.o: %.c
	$(CC) $(CXXFLAGS) $(INCLUDE) $< -c -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -c -o $@




clean:
	rm -f *.o $(TARGETS) $(TESTS) $(TESTFILES)
