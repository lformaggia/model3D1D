# ====================================================================
#   "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
#      Course on Advanced Programming for Scientific Computing
#                     Politecnico di Milano
#                         A.Y. 2014-2015
#
#                    Copyright D. Notaro 2015
# ====================================================================
#   FILE        : Makefile
#   DESCRIPTION : makefile for test simulations
#   AUTHOR      : Domenico Notaro <domenico.not@gmail.com>
#   DATE        : April 2015
# ====================================================================

include ../../../Makefile.inc

mCPPFLAGS=-I../../include -I$(mkGetFEMHome)/include
mCXXFLAGS=-std=c++11 -D=M3D1D_VERBOSE_
ifeq ($(DEBUG),yes)
  mOPTFLAGS=-g -Wall
else
  mOPTFLAGS=-O3 -march=native
  mCPPFLAGS+=-DNDEBUG
endif
mLDFLAGS= -L../../lib -L$(mkGetFEMHome)/lib
mLIBRARIES=-lproblem3d1d -lgetfem

mSRCS=$(wildcard *.cpp)
mOBJS=$(mSRCS:.cpp=.o)
mEXEC=M3D1D

mOUTDIR=vtk

.PHONY: all clean distclean

all: $(mEXEC)
	@echo
	@echo Compilation completed!

%.o: %.cpp
	$(CXX) $(mCPPFLAGS) $(mCXXFLAGS) $(mOPTFLAGS) -o $@ -c $<

$(mEXEC): $(mOBJS)
	$(CXX) $(mLDFLAGS) -o $@ $(mOBJS) $(mLIBRARIES)

clean:
	$(RM) $(mOBJS) $(mEXEC) *~ *.log

distclean: clean
	$(RM) *.txt $(mOUTDIR)/*
