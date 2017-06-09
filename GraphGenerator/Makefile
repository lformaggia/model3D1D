# ====================================================================
#   "Mixed Finite Element Methods for Coupled 3D/1D Fluid Curved Problems"
#      Course on Advanced Programming for Scientific Computing
#                     Politecnico di Milano
#                         A.Y. 2016-2017
#
#                    Copyright G. Raimondi 2017
# ====================================================================
#   FILE        : Makefile
#   DESCRIPTION : makefile for test simulations
#   AUTHOR      : Giorgio Raimondi <giorgio3.raimondi@mail.polimi.it>
#   DATE        : Maggio 2017
# ====================================================================

.PHONY: all doc pdf clean distclean library exec help


all: library
	$(MAKE) -C src/SingleBranch

help:
	@echo "make all:------ compile and create static library, and create executables"
	@echo "make help:----- print help"
	@echo "make doc:------ create doxigen file for documentation"
	@echo "make pdf:------ create doxigen file for documentation and pdf file"
	@echo "make exec:----- create in src/ the executable a run it "
	@echo "make clean:---- delete all object files"
	@echo "make distclean: delete all object files, library, executables and documentation"

library: 
	$(MAKE) -C include

doc:
	install -d doc
	doxygen Doxyfile

pdf: doc
	$(MAKE) -C doc/latex pdf
	
exec: all
	makedir vtk -C src/SingleBranch
	./src/SingleBranch/M3D1D src/SingleBranch/input.param

pdf: doc
	$(MAKE) -C doc/latex pdf

clean:
	$(RM) -r *~ *.log
	$(MAKE) clean -C include 
	$(MAKE) clean -C src/SingleBranch 

distclean: clean
	$(RM) -r doc/*
	$(MAKE) distclean -C include
	$(MAKE) distclean -C src/SingleBranch

