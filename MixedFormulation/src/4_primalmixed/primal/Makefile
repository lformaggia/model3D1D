

WARN =  -Wall -W -fmessage-length=0 -ftemplate-depth-40 -pedantic -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -Wconversion -Wredundant-decls -Wno-long-long
GETFEM_LIB =/home/nico/Downloads/getfem-4.2/src/.libs/libgetfem.a
LD_TAGS = $(GETFEM_LIB) -rdynamic /usr/lib/libblas.so.3gf /usr/lib/liblapack.so.3gf   /usr/lib/x86_64-linux-gnu/libqhull.so.6
OPT_TAGS = -O3
INCLUDE = -I/home/nico/Downloads/getfem-4.2
DEF_TAGS = -DHAVE_CONFIG -DGMM_USES_BLAS
LOG_LINKER = link.log
LOG_COMPILER = compile.log
SRC = main.cc utilities.cc import_gmshv2.cc import_network.cc
OBJ = main.o 

MCirc:	$(OBJ)
	g++ $(WARN) $(DEF_TAGS) $(OPT_TAGS) $(INCLUDE) \
	-o $@ $< $(LD_TAGS) 2>$(LOG_LINKER)

main.o:	$(SRC)
	g++ $(WARN) $(DEF_TAGS) $(OPT_TAGS) $(INCLUDE) \
	-c -o $@ $< 2>$(LOG_COMPILER)

clean:
	$(RM) *~ *.o

