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

.PHONY: all doc pdf clean distclean library exec 

all: library
	$(MAKE) -C src/SingleBranch

library: 
	$(MAKE) -C include

doc:
	install -d doc
	doxygen Doxyfile

pdf: doc
	$(MAKE) -C doc/latex pdf
	
exec: all
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

