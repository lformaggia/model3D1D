# ====================================================================
#   "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
#      Course on Advanced Programming for Scientific Computing
#                     Politecnico di Milano
#                         A.Y. 2014-2015
#
#                    Copyright D. Notaro 2015
# ====================================================================
#   FILE        : Makefile
#   DESCRIPTION : makefile for library installation
#   AUTHOR      : Domenico Notaro <domenico.not@gmail.com>
#   DATE        : February 2016
# ====================================================================

include ../../Makefile.inc

mCPPFLAGS=-I. -I$(mkGetFEMHome)/include
mCXXFLAGS=-std=c++11
mOPTFLAGS=-O3 -DNDEBUG -march=native

mLIBDIR=../lib
mLIBNAME=problem3d1d
mLIBFILE=lib$(mLIBNAME).a

mLIBSRC=$(wildcard *.cpp)
mLIBOBJS=$(mLIBSRC:.cpp=.o)
mLIBHEADERS=$(wildcard *.hpp)

.PHONY: all clean distclean library

all: library
	@echo
	@echo Library installed!

library: $(mLIBOBJS)
	install -d $(mLIBDIR)
	ar -r $(mLIBDIR)/$(mLIBFILE) $(mLIBOBJS)

%.o: %.cpp
	$(CXX) $(mCPPFLAGS) $(mCXXFLAGS) $(mOPTFLAGS) -o $@ -c $<

clean:
	$(RM) $(mLIBOBJS) *~

distclean: clean
	$(RM) $(mLIBDIR)/$(mLIBFILE)
