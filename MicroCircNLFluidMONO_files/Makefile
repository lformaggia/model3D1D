
GETFEM_DIR=/opt/getfem4.2
GETFEM_INCLUDE_DIR=$(GETFEM_DIR)/include
GETFEM_LIB_DIR=$(GETFEM_DIR)/lib
BOOST_DIR=/usr/local/boost_1_57_0
BOOST_INCLUDE_DIR=$(BOOST_DIR)/include
OPTFLAGS=-g
BLAS_LIB=-lblas
GETFEM_LIB=-L$(GETFEM_LIB_DIR) -lgetfem
WARNFLAGS= -Wall -W -fmessage-length=0 -ftemplate-depth-40 -pedantic -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -Wconversion -Wredundant-decls -Wno-long-long

LDFLAGS= $(OPTFLAGS) -rdynamic -Wl,-rpath,$(GETFEM_LIB_DIR)
LDLIBS= $(GETFEM_LIB) $(BLAS_LIB)
CXXFLAGS=$(OPTFLAGS) $(WARNFLAGS) $(INCLUDE)
CFLAGS=$(OPTFLAGS) $(WARNFLAGS) $(INCLUDE)
CPPFLAGS=-DHAVE_CONFIG -DGMM_USES_BLAS
INCLUDE = -I. -I$(GETFEM_INCLUDE_DIR) -I$(GETFEM_INCLUDE_DIR)/gmm   -I$(BOOST_INCLUDE_DIR)
CXX=g++-4.9
CC=$(CXX)
DEF_TAGS = 
LOG_LINKER = link.log
LOG_COMPILER = compile.log
SRC = MCirc.cc 
OBJ = MCirc.o 

MCirc: $(OBJ)

clean:
	$(RM) *~ *.o MCirc

