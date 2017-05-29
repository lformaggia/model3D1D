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

.PHONY: all doc clean distclean library exec1 

all: library


library: 
	$(MAKE) -C include
	$(MAKE) -C src/SingleBranch

doc:
	install -d doc
	doxygen Doxyfile
	
exec1: all
	makedir vtk -C src/SingleBranch
	./M3D1D input.param

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

