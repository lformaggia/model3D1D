
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
include ../../Makefile.inc

VERBOSE= -D=M3D1D_VERBOSE_

SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)
LIBHEADERS=$(wildcard *.hpp)
LIBNAME= libc_problem3d1d.a
LIBDIR = ../lib

.PHONY: all clean distclean 

all: $(LIBDIR)/$(LIBNAME)
	@echo
	@echo Compilation completed!

%.o: %.cpp
	$(CXX) $(STDFLAG) $(VERBOSE) $(C_MANWORKS_INCLUDE) $(OPTFLAGS) -o $@ -c $<

$(LIBDIR)/$(LIBNAME): $(OBJS)
	install -d $(LIBDIR)
	ar -r $@ $(OBJS)

clean:
	$(RM) $(OBJS) $(EXEC) *~ *.log

distclean: clean
	$(RM) *.txt
	$(RM) $(LIBDIR)/$(LIBNAME) 
