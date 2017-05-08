# Inlcuding general Makefile.in
include ../../../Makefile.inc
# Including specific Graph generator Makefile.inc
include ../../Makefile.inc

.PHONY= help all exec clean distclean

.DEFAUL_GOAL: all

all: exec

help:
	@echo "make help: ----- Prints this help"
	@echo "make all: ------ "
	@echo "make clean: ---- Cleans all object files"

main.o: main.cpp
	$(CXX) $(STDFLAG) $(GRAPH_INCLUDE) $(OPTFLAGS) $(WARNREAD)  -o main.o -c main.cpp

exec: main.o
	$(CXX) $(STDFLAG) $(GRAPH_LIBRARY) -L../../lib -lgraphgenerator  main.o -o main  

clean:
	$(RM) $(LIBOBJS) *~

distclea:
	$(RM) main
