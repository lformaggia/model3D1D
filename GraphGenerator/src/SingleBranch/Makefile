# Inlcuding general Makefile.in
include ../../../Makefile.inc
# Including specific Graph generator Makefile.inc
include ../../Makefile.inc

LDFLAGS= -L/home/pacs_student/model3D1D/GraphGenerator/lib  $(GRAPH_LIBRARY_DIR)
LIB= -lgraphgenerator $(GRAPH_LIBRARY)
EXEC=main
.PHONY= help all clean distclean run

.DEFAUL_GOAL: all

all: $(EXEC)

help:
	@echo "make help: ----- Prints this help"
	@echo "make all: ------ Complile all object and link libraries"
	@echo "make clean: ---- Cleans all object files"
	@echo "make distclean:- Cleans all object and executables"

main.o: main.cpp
	$(CXX) $(STDFLAG) $(GRAPH_INCLUDE) $(OPTFLAGS) $(WARNREAD)  -o main.o -c main.cpp

 $(EXEC): main.o
	$(CXX)  $(LDFLAGS)  -o $(EXEC) main.o $(LIB)

run:  $(EXEC)
	install -d vtk
	./main input.param
clean:
	$(RM) $(LIBOBJS) *~

distclean: clean
	$(RM) main *.pts
