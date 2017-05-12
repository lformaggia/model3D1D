# Inlcuding general Makefile.in
include ../../../Makefile.inc
# Including specific Graph generator Makefile.inc
include ../../Makefile.inc

LIB=-L/home/pacs_student/model3D1D/GraphGenerator/lib   -lgraphgenerator $(GRAPH_LIBRARY)
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
	$(CXX)    -o main main.o $(LIB)

clean:
	$(RM) $(LIBOBJS) *~

distclea:
	$(RM) main
