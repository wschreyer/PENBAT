###################################################################
# This Makefile was created using the bat-project script
# for project PENBAT
# bat-project is part of Bayesian Analysis Toolkit (BAT).
# BAT can be downloaded from http://mpp.mpg.de/bat
###################################################################
#
# Run 'make' to compile the program and 'make clean' to remove
# all compiled parts and 'clean' the directory.
#
# You might need to adjust the CXXFLAGS and LIBS based on
# the BAT installation on your system. Consult the gmake manual
# for details.
#
###################################################################

# compiler and flags
CXX          = g++
CXXFLAGS     = -O3 -Wall -fPIC -Wno-deprecated
LD           = /usr/bin/ld -m elf_x86_64
LDFLAGS      = -O3

# ----------------------------------------------------------------------
# The following definitions rely on the script bat-config being
# available in $PATH. If BAT is not installed in the standard system
# directories, update $PATH accordingly.

CXXFLAGS += `bat-config --cflags`
LIBS := `bat-config --libs`

# List of all classes (models) used in the program
# Add classes to the end. Backslash indicates continuation
# on the next line
CXXSRCS      = \
	runPENModel.cxx PENModel.cxx

# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#

CXXOBJS      = $(patsubst %.cxx,%.o,$(CXXSRCS))
MYPROGS     = runPENBAT

GARBAGE      = $(CXXOBJS) *.o *~ link.d $(MYPROGS)

# targets
all : runPENModel

link.d : $(patsubst %.cxx,%.h,$(CXXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d;

-include link.d

%.o : %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	rm -f $(GARBAGE)

runPENModel : $(CXXOBJS)
	$(CXX) $(LDFLAGS) $(CXXOBJS) $(LIBS) -o runPENBAT

print :
	@echo compiler  : $(CXX)
	@echo c++ srcs  : $(CXXSRCS)
	@echo c++ objs  : $(CXXOBJS)
	@echo c++ flags : $(CXXFLAGS)
	@echo ld flags  : $(LDFLAGS)
	@echo libs      : $(LIBS)

	
analyze : analyze.cxx
	$(CXX) `root-config --cflags` `root-config --libs` analyze.cxx -o analyze
	
simulate : simulate.cxx
	$(CXX) -std=c++11 `bat-config --cflags` `bat-config --libs` simulate.cxx -o simulate
