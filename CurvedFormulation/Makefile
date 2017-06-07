#======================================================================
#    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Curved Problems"
#        Course on Advanced Programming for Scientific Computing
#                      Politecnico di Milano
#                          A.Y. 2016-2017
#                  
#               Copyright (C) 2017 Giorgio Raimondi
#======================================================================
# 
#  @file   Makefile
#  @author Giorgio Raimondi <giorgio3.raimondi@mail.polimi.it>
#  @date   May 2017.
#
#======================================================================    
.PHONY= all clean distclean doc pdf library
.DEFAULT_GOAL=all



all:
	make -C include
	make -C src/1_LinearSingleBranch
	make -C src/2_CurvedSingleBranch
	make -C src/3_LinearBifurcation
	make -C src/4b_CurvedBifurcation
	make -C src/4_CurvedBifurcation
	make -C src/5_LinearRombus
	make -C src/6_CurvedRombus
	make -C src/6b_Circular_Rombus


library:
	make -C include

doc:
	install -d doc
	doxygen Doxyfile

pdf: doc
	$(MAKE) -C doc/latex pdf

clean:
	make clean -C include
	make clean -C src/1_LinearSingleBranch
	make clean -C src/2_CurvedSingleBranch
	make clean -C src/3_LinearBifurcation
	make clean -C src/4_CurvedBifurcation
	make clean -C src/4b_CurvedBifurcation
	make clean -C src/5_LinearRombus
	make clean -C src/6_CurvedRombus
	make clean -C src/6b_Circular_Rombus

distclean: clean
	$(RM) -r doc/*
	$(RM) -r doc
	make distclean -C include
	make distclean -C src/1_LinearSingleBranch
	make distclean -C src/2_CurvedSingleBranch
	make distclean -C src/3_LinearBifurcation
	make distclean -C src/4_CurvedBifurcation
	make distclean -C src/4b_CurvedBifurcation
	make distclean -C src/5_LinearRombus
	make distclean -C src/6_CurvedRombus
	make distclean -C src/6b_Circular_Rombus
