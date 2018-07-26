PROGRAM	= Kprocessor

CC	= gcc-6
CXX	= g++-6



INCLUDE	= -IThirdParty/seqan/include -IMQF
CPPFLAGS	= -Wall -Wextra -std=c++14 -fPIC -fopenmp -W -Wall -pedantic
LDFLAGS	= -lrt -lpthread -lbz2 -lz


SEQAN_FLAGS= -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1  -DSEQAN_HAS_OPENMP=1 -Wl,--whole-archive -Wl,--no-whole-archive
CPPFLAGS += $(SEQAN_FLAGS)
# The directories in which source files reside.
# If not specified, all subdirectories of the current directory will be added recursively.
SRCDIRS	:=  KmerCounter/ HashUtils/
UNAME_S  := $(shell uname -s)


#Actually $(INCLUDE) is included in $(CPPFLAGS).
CPPFLAGS      += $(INCLUDE)

## Implicit Section: change the following only when necessary.
##==========================================================================
# The source file types (headers excluded).
# .c indicates C source files, and others C++ ones.
SRCEXTS = .c .C .cc .cpp .CPP .c++ .cxx .cp

# The header file types.
HDREXTS = .h .H .hh .hpp .HPP .h++ .hxx .hp

# The pre-processor and compiler options.
# Users can override those variables from the command line.
CFLAGS  = -O3
CXXFLAGS= -O3

# The command used to delete file.
RM     = rm -f

ETAGS = etags
ETAGSFLAGS =

CTAGS = ctags
CTAGSFLAGS =

## Stable Section: usually no need to be changed. But you can add more.
##==========================================================================
ifeq ($(SRCDIRS),)
	SRCDIRS := $(shell find $(SRCDIRS) -type d)
endif
SOURCES = $(foreach d,$(SRCDIRS),$(wildcard $(addprefix $(d)/*,$(SRCEXTS))))
HEADERS = $(foreach d,$(SRCDIRS),$(wildcard $(addprefix $(d)/*,$(HDREXTS))))
SRC_CXX = $(filter-out %.c,$(SOURCES))
OBJS	= counterMain.o estimateMemoryMain.o dumpMain.o
OBJS    += $(addsuffix .o, $(basename $(SOURCES)))
OBJS += MQF/gqf.o MQF/utils.o
#DEPS    = $(OBJS:%.o=%.d) #replace %.d with .%.d (hide dependency files)
DEPS    = $(foreach f, $(OBJS), $(addprefix $(dir $(f))., $(patsubst %.o, %.d, $(notdir $(f)))))


## Define some useful variables.
DEP_OPT = $(shell if `$(CC) --version | grep -i "GCC" >/dev/null`; then \
                  echo "-MM"; else echo "-M"; fi )
DEPEND.d    = $(CC)  $(DEP_OPT)  $(EXTRA_CFLAGS) $(CFLAGS) $(CPPFLAGS)
COMPILE.c   = $(CC)  $(EXTRA_CFLAGS) $(CFLAGS)   $(CPPFLAGS) -c
COMPILE.cxx = $(CXX) $(EXTRA_CFLAGS) $(CXXFLAGS) $(CPPFLAGS) -c
LINK.c      = $(CC)  $(EXTRA_CFLAGS) $(CFLAGS)   $(CPPFLAGS) $(LDFLAGS)
LINK.cxx    = $(CXX) $(EXTRA_CFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS)

.PHONY: all objs tags ctags clean distclean help show

# Delete the default suffixes
.SUFFIXES:



all: $(PROGRAM)

MQF/gqf.o:
	cd MQF && make libgqf
MQF/utils.o:
		cd MQF && make libgqf
# Rules for creating dependency files (.d).
#------------------------------------------

.%.d:%.c
	@echo -n $(dir $<) > $@
	@$(DEPEND.d) $< >> $@

.%.d:%.C
	@echo -n $(dir $<) > $@
	@$(DEPEND.d) $< >> $@

.%.d:%.cc
	@echo -n $(dir $<) > $@
	@$(DEPEND.d) $< >> $@

.%.d:%.cpp
	@echo -n $(dir $<) > $@
	@$(DEPEND.d) $< >> $@

.%.d:%.CPP
	@echo -n $(dir $<) > $@
	@$(DEPEND.d) $< >> $@

.%.d:%.c++
	@echo -n $(dir $<) > $@
	@$(DEPEND.d) $< >> $@

.%.d:%.cp
	@echo -n $(dir $<) > $@
	@$(DEPEND.d) $< >> $@

.%.d:%.cxx
	@echo -n $(dir $<) > $@
	@$(DEPEND.d) $< >> $@

# Rules for generating object files (.o).
#----------------------------------------
objs:$(OBJS)

%.o:%.c
	$(COMPILE.c) $< -o $@

%.o:%.C
	$(COMPILE.cxx) $< -o $@

%.o:%.cc
	$(COMPILE.cxx) $< -o $@

%.o:%.cpp
	$(COMPILE.cxx) $< -o $@

%.o:%.CPP
	$(COMPILE.cxx) $< -o $@

%.o:%.c++
	$(COMPILE.cxx) $< -o $@

%.o:%.cp
	$(COMPILE.cxx) $< -o $@

%.o:%.cxx
	$(COMPILE.cxx) $< -o $@

# Rules for generating the tags.
#-------------------------------------
tags: $(HEADERS) $(SOURCES)
	$(ETAGS) $(ETAGSFLAGS) $(HEADERS) $(SOURCES)

ctags: $(HEADERS) $(SOURCES)
	$(CTAGS) $(CTAGSFLAGS) $(HEADERS) $(SOURCES)

# Rules for generating the executable.
#-------------------------------------
$(PROGRAM):$(OBJS) main.o
	$(LINK.cxx)  main.o $(OBJS) $(EXTRA_LDFLAGS) $(LDFLAGS) -o $@
	@echo Type ./$@ to execute the program.

ifndef NODEP
ifneq ($(DEPS),)
  sinclude $(DEPS)
endif
endif


TESTS = tests/testsMain.o tests/testKmerCounter.o
test: $(OBJS) $(TESTS)
	$(LINK.cxx) $(EXTRA_LDFLAGS)  $^ $(LDFLAGS) -o $@
	rm -f tests/testData/tmp*


clean:
	$(RM) $(OBJS) $(PROGRAM) $(PROGRAM).exe $(TESTS)
	rm -f tests/testData/tmp*
	cd MQF && make clean

distclean: clean
	$(RM) $(DEPS) TAGS





# Show help.
help:
	@echo "Kprocessor Makefile"
	@echo
	@echo 'Usage: make [TARGET]'
	@echo 'TARGETS:'
	@echo '  all       (=make) compile and link.'
	@echo '  NODEP=yes make without generating dependencies.'
	@echo '  objs      compile only (no linking).'
	@echo '  tags      create tags for Emacs editor.'
	@echo '  ctags     create ctags for VI editor.'
	@echo '  clean     clean objects and the executable file.'
	@echo '  distclean clean objects, the executable and dependencies.'
	@echo '  show      show variables (for debug use only).'
	@echo '  help      print this message.'
	@echo

# Show variables (for debug use only.)
show:
	@echo 'PROGRAM     :' $(PROGRAM)
	@echo 'SRCDIRS     :' $(SRCDIRS)
	@echo 'HEADERS     :' $(HEADERS)
	@echo 'SOURCES     :' $(SOURCES)
	@echo 'SRC_CXX     :' $(SRC_CXX)
	@echo 'OBJS        :' $(OBJS)
	@echo 'DEPS        :' $(DEPS)
	@echo 'DEPEND      :' $(DEPEND)
	@echo 'DEPEND.d    :' $(DEPEND.d)
	@echo 'COMPILE.c   :' $(COMPILE.c)
	@echo 'COMPILE.cxx :' $(COMPILE.cxx)
	@echo 'link.c      :' $(LINK.c)
	@echo 'link.cxx    :' $(LINK.cxx)

## End of the Makefile ##  Suggestions are welcome  ## All rights reserved ##
#############################################################################

## Thanks Pear, kevin1078 , and whyglinux for the Makefile template ##
## https://github.com/Cheedoong/MakefileTemplate                    ##
######################################################################
