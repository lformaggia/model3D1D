

DEFAULT_GOAL= all

.PHONY= all

all:
	make -C MixedFormulation
	make -C CurvedFormulation
	make -C GraphGenerator

clean: 
	make clean -C MixedFormulation
	make clean -C CurvedFormulation
	make clean -C GraphGenerator

distclean: 
	make distclean -C MixedFormulation
	make distclean -C CurvedFormulation
	make distclean -C GraphGenerator
