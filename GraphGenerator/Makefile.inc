CPPFLAGS= $(GRAPH_INCLUDES)
OPTFLAGS=-O3 -DNDEBUG -march=native

LIBDIR=../lib
LIBNAME=graphgenerator
LIBFILE_STATIC=lib$(LIBNAME).a

LIBSRC=$(wildcard *.cpp)
LIBOBJS=$(LIBSRC:.cpp=.o)
WARNREAD=$(WARN) -Wno-deprecated -Wunknown-pragmas
