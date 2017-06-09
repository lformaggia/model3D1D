# Including generic Makefile.inc
include ../../../Makefile.inc

LIB= -lc_problem3d1d $(C_MANWORKS_LIBRARY)  
LIBDIR= $(C_MANWORK_LIBRARY_DIR) -L$(C_MANWORKS_LIB_DIR) 
.PHONY= help all main run clean distclean

.DEFAUL_GOAL: all

all: main
	@echo 
	@echo "3_LinearBifurcation test case compiled"
	@echo

help:
	@echo "make help: ------- Prints this help"
	@echo "make all: -------- Compile and prints help"
	@echo "make main: ------- Compile the program"
	@echo "make run: -------- Run the executable with input.param file"
	@echo "make clean: ------ Cleans all object files"
	@echo "make distclean: -- Cleans all object files and remove executable"

main.o: main.cpp
	$(CXX) $(STDFLAG) $(C_MANWORKS_INCLUDE) $(OPTFLAGS) $(WARN)  -o main.o -c main.cpp

main: main.o
	$(CXX) $(LIBDIR) -o main main.o  $(LIB)

run: main
	./main input.param

clean:
	$(RM) main.o *~

distclean:
	$(RM) main



