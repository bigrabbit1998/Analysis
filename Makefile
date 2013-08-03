# Need this to get SHAREDSUFFIX 
-include ../config.mk

# A few variables used in this Makefile:
EX           := hist tree
EXE          := $(addsuffix .exe,$(EX))
PYTHIA8      ?= $(PWD)/..
STATICLIB    := $(PYTHIA8)/lib/archive/libpythia8.a
SHAREDLIB    := $(PYTHIA8)/lib/libpythia8.$(SHAREDSUFFIX)
DICTCXXFLAGS := -I$(PYTHIA8)/include 
ROOTCXXFLAGS := $(DICTCXXFLAGS) $(shell root-config --cflags)

# Libraries to include 
ifeq (x$(ENABLEGZIP),xyes)
LIBGZIP=-L$(BOOSTLIBLOCATION) -lboost_iostreams -L$(ZLIBLOCATION) -lz
endif

# LDFLAGS1 for static library, LDFLAGS2 for shared library
LDFLAGS1 := $(shell root-config --ldflags --glibs) \
  -L$(PYTHIA8)/lib/archive -lpythia8 -llhapdfdummy $(LIBGZIP)
LDFLAGS2 := $(shell root-config --ldflags --glibs) \
  -L$(PYTHIA8)/lib -lpythia8 -llhapdfdummy $(LIBGZIP)

all: $(EX)

# fastjet
FASTJET     = $(shell fastjet-config --prefix)
FASTJETLIBS = -L$(FASTJET)/lib -lfastjet
FASTJETINCS = $(shell fastjet-config --cxxflags)

CXXFLAGS   += $(FASTJETINCS) $(ROOTCXXFLAGS)
OUF        += fjClustering.o MyEvent.o JetSeparation.o MyTopEvent.o 


install : all
all : hist
# building rules
hist: $(STATICLIB) hist.cc $(OUF)
	$(CXX) $(CXXFLAGS) $@.cc -g $(OUF) -o $@.exe $(LDFLAGS1) $(FASTJETLIBS)

mydat: $(STATICLIB) mydat.cxx $(OUF)
	$(CXX) $(CXXFLAGS) $@.cxx -g $(OUF) -o $@.exe $(LDFLAGS1) $(FASTJETLIBS)

# Rule to build tree 
tree: $(STATICLIB) tree.cc
	rootcint -f treeDict.cc -c $(DICTCXXFLAGS) pythiaROOT.h pythiaLinkdef.h
	$(CXX) $(ROOTCXXFLAGS) treeDict.cc $@.cc -g -o $@.exe $(LDFLAGS1)

.cxx.o:  
	$(CXX) $(CXXFLAGS) -c $< 
.cpp.o:  
	$(CXX) $(CXXFLAGS) -c $< 


# Rule to build full dictionary
dict: $(SHAREDLIB)
	rootcint -f pythiaDict.cc -c $(DICTCXXFLAGS) \
           -DPYTHIA8_COMPLETE_ROOT_DICTIONARY \
           pythiaROOT.h pythiaLinkdef.h
	$(CXX) -shared -fPIC -g -o pythiaDict.$(SHAREDSUFFIX) pythiaDict.cc \
         -DPYTHIA8_COMPLETE_ROOT_DICTIONARY \
         $(ROOTCXXFLAGS) $(LDFLAGS2)


# Error messages if PYTHIA libraries don't exist
$(STATICLIB):
	@echo "Error: PYTHIA 8 archive library must be built"
	@false
$(SHAREDLIB):
	@echo "Error: PYTHIA 8 shared library must be built"
	@false

# Clean up
clean:
	rm -rf ./.libs ./.obj *.o *.la *~ $(EXE) hist.root pythiaDict.* \
		reeDict.cc treeDict.h pytree.root

