# Inlcuding general Makefile.in
include ../../Makefile.inc
# Including specific Graph generator Makefile.inc
include ../Makefile.inc

.PHONY= help all exec clean distclean

.DEFAUL_GOAL: all


all: static_lib
	@echo
	@echo "The libraries installation is completed!!!"
	@echo

help:
	@echo "make help: ----- Prints this help"
	@echo "make all: ------ Makes static and dynamic library"
	@echo "make clean: ---- Cleans all object files"
	@echo "make distclean:  Cleans all"

static_lib: $(LIBOBJS)
	install -d $(LIBDIR)
	ar -r $(LIBDIR)/$(LIBFILE_STATIC) $(LIBOBJS)

%.o: %.cpp
	$(CXX) $(STDFLAG) $(CPPFLAGS) $(GRAPH_INCLUDE) $(OPTFLAGS) $(WARNREAD) -o $@ -c $<

clean:
	$(RM) $(LIBOBJS) *~

distclean: clean
	$(RM) $(LIBDIR)/$(LIBFILE_STATIC)
